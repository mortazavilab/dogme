#!/usr/bin/env python3
"""Add Splitcode cell-barcode and UMI tags to an annotated BAM by read ID."""

from __future__ import annotations

import argparse
import gzip
import sqlite3
import tempfile
from pathlib import Path
from typing import Iterator, TextIO


def _open_fastq(path: Path) -> TextIO:
    if path.name.endswith((".gz", ".bgz", ".bgzip")):
        return gzip.open(path, "rt")
    return path.open()


def _read_fastq(path: Path) -> Iterator[tuple[str, str]]:
    with _open_fastq(path) as fastq:
        line_number = 0
        while True:
            header = fastq.readline()
            if not header:
                return
            sequence = fastq.readline()
            separator = fastq.readline()
            quality = fastq.readline()
            line_number += 4
            if not sequence or not separator or not quality or not header.startswith("@") or not separator.startswith("+"):
                raise ValueError(f"Invalid FASTQ record in {path} ending at line {line_number}")
            sequence = sequence.rstrip("\r\n")
            quality = quality.rstrip("\r\n")
            if not sequence or len(sequence) != len(quality):
                raise ValueError(
                    f"Empty sequence or sequence-quality length mismatch in {path} at line {line_number - 2}"
                )
            read_id = header[1:].split()[0]
            if not read_id:
                raise ValueError(f"Empty FASTQ read ID in {path} at line {line_number - 3}")
            yield read_id, sequence.upper()


def tag_bam(
    bam_path: Path,
    umi_path: Path,
    barcode_path: Path,
    output_path: Path,
    qc_output_path: Path,
) -> dict[str, int]:
    try:
        import pysam
    except ImportError as exc:
        raise RuntimeError("BAM tagging requires pysam in the DOGME execution image.") from exc

    for path in (bam_path, umi_path, barcode_path):
        if not path.is_file():
            raise FileNotFoundError(f"Required tagging input does not exist: {path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    qc_output_path.parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="dogme_bam_tags_") as temporary_directory:
        db = sqlite3.connect(Path(temporary_directory) / "read_tags.sqlite")
        try:
            db.execute(
                "CREATE TABLE read_tags "
                "(read_id TEXT PRIMARY KEY, umi TEXT NOT NULL, barcode TEXT NOT NULL, matched INTEGER NOT NULL DEFAULT 0)"
            )
            umi_records = iter(_read_fastq(umi_path))
            barcode_records = iter(_read_fastq(barcode_path))
            fastq_count = 0
            while True:
                umi_record = next(umi_records, None)
                barcode_record = next(barcode_records, None)
                if umi_record is None and barcode_record is None:
                    break
                if umi_record is None or barcode_record is None:
                    raise ValueError("Splitcode UMI and barcode FASTQs contain different record counts")
                umi_id, umi = umi_record
                barcode_id, barcode = barcode_record
                if umi_id != barcode_id:
                    raise ValueError(
                        f"Splitcode FASTQ read IDs are out of sync: {umi_id!r} != {barcode_id!r}"
                    )
                try:
                    db.execute(
                        "INSERT INTO read_tags (read_id, umi, barcode) VALUES (?, ?, ?)",
                        (umi_id, umi, barcode),
                    )
                except sqlite3.IntegrityError as exc:
                    raise ValueError(f"Duplicate read ID in Splitcode FASTQs: {umi_id}") from exc
                fastq_count += 1
            db.commit()
            if not fastq_count:
                raise ValueError("Splitcode FASTQs contain no records")

            total_alignments = 0
            tagged_alignments = 0
            matched_read_ids = 0
            with pysam.AlignmentFile(str(bam_path), "rb") as source:
                with pysam.AlignmentFile(str(output_path), "wb", template=source) as target:
                    for read in source.fetch(until_eof=True):
                        total_alignments += 1
                        tags = db.execute(
                            "SELECT umi, barcode FROM read_tags WHERE read_id = ?",
                            (read.query_name,),
                        ).fetchone()
                        if tags is not None:
                            umi, barcode = tags
                            read.set_tag("CB", barcode, value_type="Z", replace=True)
                            read.set_tag("UB", umi, value_type="Z", replace=True)
                            db.execute(
                                "UPDATE read_tags SET matched = 1 WHERE read_id = ?",
                                (read.query_name,),
                            )
                            tagged_alignments += 1
                        target.write(read)
            db.commit()
            missing_ids = db.execute(
                "SELECT read_id FROM read_tags WHERE matched = 0 ORDER BY read_id LIMIT 10"
            ).fetchall()
            missing_count = db.execute(
                "SELECT COUNT(*) FROM read_tags WHERE matched = 0"
            ).fetchone()[0]
            if missing_count:
                examples = ", ".join(read_id for (read_id,) in missing_ids)
                raise ValueError(
                    f"{missing_count} Splitcode read IDs were absent from the annotated BAM"
                    f" (examples: {examples})"
                )
        finally:
            db.close()

    metrics = {
        "splitcode_reads": fastq_count,
        "bam_alignments": total_alignments,
        "tagged_bam_alignments": tagged_alignments,
        "untagged_bam_alignments": total_alignments - tagged_alignments,
    }
    with qc_output_path.open("w") as output:
        output.write("metric\tvalue\n")
        for metric, value in metrics.items():
            output.write(f"{metric}\t{value}\n")
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path)
    parser.add_argument("--umi", required=True, type=Path)
    parser.add_argument("--barcode", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--qc-output", required=True, type=Path)
    args = parser.parse_args()
    tag_bam(args.bam, args.umi, args.barcode, args.output, args.qc_output)


if __name__ == "__main__":
    main()
