#!/usr/bin/env python3
"""Normalize and combine splitcode FASTQs for single-cell quantification."""

import argparse
import gzip
from pathlib import Path
from typing import TextIO


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


def write_record(output: TextIO, record: list[str], orientation: str, fastq_type: str) -> None:
    if len(record) != 4 or not record[0].startswith("@") or not record[2].startswith("+"):
        raise ValueError("invalid FASTQ record")
    if fastq_type == "barcode":
        record[1] = transform_barcode(record[1], orientation) + "\n"
    output.writelines(record)


def combine_fastqs(sample: str) -> None:
    for fastq_type in FASTQ_TYPES:
        output_path = Path(f"{sample}_{fastq_type}.fastq.gz")
        with open_fastq(output_path, "w") as output:
            for orientation in ORIENTATIONS:
                input_path = Path(f"{orientation}_{fastq_type}.fastq.gz")
                with open_fastq(input_path, "r") as input_file:
                    while record := [input_file.readline() for _ in range(4)]:
                        if not record[0]:
                            break
                        write_record(output, record, orientation, fastq_type)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", required=True)
    args = parser.parse_args()
    combine_fastqs(args.sample)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())