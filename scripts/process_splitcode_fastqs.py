#!/usr/bin/env python3
"""Normalize and combine splitcode FASTQs for single-cell quantification."""

import argparse
import gzip
from pathlib import Path
from typing import Iterator, TextIO


FASTQ_TYPES = ("barcode", "cDNA", "umi")
ORIENTATIONS = ("f", "c", "r", "rc")


def open_fastq(path: Path, mode: str) -> TextIO:
    return gzip.open(path, mode + "t")


def transform_barcode(sequence: str, orientation: str) -> str:
    sequence = sequence.rstrip("\r\n")
    if orientation == "c":
        return "".join(sequence[start : start + 8][::-1] for start in range(0, 24, 8))
    if orientation == "r":
        return sequence[::-1]
    if orientation == "rc":
        return sequence[16:24] + sequence[8:16] + sequence[0:8]
    return sequence


def transform_quality(quality: str, orientation: str) -> str:
    if orientation == "c":
        return "".join(quality[start : start + 8][::-1] for start in range(0, 24, 8))
    if orientation == "r":
        return quality[::-1]
    if orientation == "rc":
        return quality[16:24] + quality[8:16] + quality[0:8]
    return quality


def write_record(output: TextIO, record: list[str], orientation: str, fastq_type: str) -> None:
    if len(record) != 4 or not record[0].startswith("@") or not record[2].startswith("+"):
        raise ValueError("invalid FASTQ record")
    if fastq_type == "barcode":
        record[1] = transform_barcode(record[1], orientation) + "\n"
        record[3] = transform_quality(record[3].rstrip("\r\n"), orientation) + "\n"
    output.writelines(record)


def read_records(path: Path) -> Iterator[list[str]]:
    with open_fastq(path, "r") as input_file:
        while True:
            record = [input_file.readline() for _ in range(4)]
            if not record[0]:
                return
            if len(record) != 4 or not record[3] or not record[2].startswith("+"):
                raise ValueError(f"invalid FASTQ record in {path}")
            if len(record[1].rstrip("\r\n")) != len(record[3].rstrip("\r\n")):
                raise ValueError(f"sequence-quality length mismatch in {path}")
            yield record


def orientation_paths(orientation: str) -> dict[str, Path]:
    return {
        fastq_type: Path(f"{orientation}_{fastq_type}.fastq.gz")
        for fastq_type in FASTQ_TYPES
    }


def read_triplets(orientation: str) -> Iterator[dict[str, list[str]]]:
    paths = orientation_paths(orientation)
    if not all(path.is_file() for path in paths.values()):
        return
    inputs = {
        fastq_type: read_records(path) for fastq_type, path in paths.items()
    }
    seen_ids: set[str] = set()
    while True:
        records = {
            fastq_type: next(iterator, None)
            for fastq_type, iterator in inputs.items()
        }
        if all(record is None for record in records.values()):
            return
        if any(record is None for record in records.values()):
            raise ValueError(f"orientation {orientation} has unsynchronized FASTQs")
        identifiers = {
            fastq_type: records[fastq_type][0].split()[0]
            for fastq_type in FASTQ_TYPES
        }
        if len(set(identifiers.values())) != 1:
            raise ValueError(f"orientation {orientation} has mismatched read identifiers")
        read_id = next(iter(identifiers.values()))
        if read_id in seen_ids:
            raise ValueError(f"duplicate read identifier within orientation: {read_id}")
        seen_ids.add(read_id)
        yield records


def write_qc(path: Path, values: dict[str, str | int]) -> None:
    with path.open("w") as output:
        output.write("metric\tvalue\n")
        for key, value in values.items():
            output.write(f"{key}\t{value}\n")


def combine_fastqs(sample: str, qc_output: Path | None = None) -> None:
    output_paths = {
        fastq_type: Path(f"{sample}_{fastq_type}.fastq.gz") for fastq_type in FASTQ_TYPES
    }
    outputs = {fastq_type: open_fastq(path, "w") for fastq_type, path in output_paths.items()}
    present_orientations = [
        orientation
        for orientation in ORIENTATIONS
        if all(path.is_file() for path in orientation_paths(orientation).values())
    ]
    if not present_orientations:
        raise FileNotFoundError("no complete splitcode orientation FASTQ set found")
    id_orientations: dict[str, set[str]] = {}
    empty_records = 0
    input_records = 0
    try:
        for orientation in present_orientations:
            for records in read_triplets(orientation):
                input_records += 1
                read_id = records["cDNA"][0].split()[0]
                if any(not records[fastq_type][1].strip() for fastq_type in FASTQ_TYPES):
                    empty_records += 1
                    continue
                id_orientations.setdefault(read_id, set()).add(orientation)

        ambiguous_ids = {
            read_id for read_id, orientations in id_orientations.items() if len(orientations) > 1
        }
        emitted_records = 0
        for orientation in present_orientations:
            for records in read_triplets(orientation):
                read_id = records["cDNA"][0].split()[0]
                if read_id in ambiguous_ids or any(
                    not records[fastq_type][1].strip() for fastq_type in FASTQ_TYPES
                ):
                    continue
                for fastq_type in FASTQ_TYPES:
                    write_record(outputs[fastq_type], records[fastq_type], orientation, fastq_type)
                emitted_records += 1
        if qc_output is not None:
            write_qc(qc_output, {
                "orientations_present": ",".join(present_orientations),
                "orientations_missing": ",".join(
                    orientation for orientation in ORIENTATIONS if orientation not in present_orientations
                ),
                "input_triplets": input_records,
                "empty_triplets_excluded": empty_records,
                "ambiguous_read_ids_excluded": len(ambiguous_ids),
                "emitted_triplets": emitted_records,
            })
    finally:
        for output in outputs.values():
            output.close()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--qc-output", type=Path)
    args = parser.parse_args()
    combine_fastqs(args.sample, args.qc_output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())