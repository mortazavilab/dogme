#!/usr/bin/env python3
import hashlib
import os
import argparse
import csv
from pathlib import Path
from datetime import datetime

def compute_md5(file_path):
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def collect_file_info(directory, extensions):
    file_data = []
    for subfolder, ext in [("bams", ".bam"), ("bedMethyl", ".bed"), ("openChromatin", ".bed")]:
        target_dir = Path(directory) / subfolder
        if not target_dir.exists():
            continue
        for filepath in target_dir.rglob(f'*{ext}'):
            size_bytes = filepath.stat().st_size
            file_data.append({
                'filename': filepath.name,
                'extension': ext,
                'size_bytes': size_bytes,
                'path': str(filepath.resolve()),
                'mod_time': datetime.fromtimestamp(filepath.stat().st_mtime),
                'md5sum': compute_md5(filepath)
            })
    return file_data

def write_report(data, output_file):
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['filename', 'extension', 'size_bytes', 'path', 'mod_time', 'md5sum']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(data)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True, help="Directory to scan for output files")
    parser.add_argument("-o", "--output_report", required=True, help="Output report file path (.tsv)")
    args = parser.parse_args()

    extensions = ['.bam', '.bed']
    data = collect_file_info(args.input_dir, extensions)
    write_report(data, args.output_report)

if __name__ == "__main__":
    main()
