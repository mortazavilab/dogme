#!/usr/bin/env python3
# This script generates two reports: a file inventory (TSV) and a QC summary (CSV).

import argparse
import csv
import glob
import gzip
import hashlib
import os
import sys
from datetime import datetime
from pathlib import Path
from statistics import median

try:
    import pysam
except ImportError:
    pysam = None

# ==============================================================================
# SECTION 1: FILE INVENTORY FUNCTIONS
# ==============================================================================

def compute_md5(file_path):
    """Calculates the MD5 checksum of a file."""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()

def collect_file_info(directory):
    """Collects file metadata (size, mod_time, md5) from subdirectories."""
    file_data = []
    search_patterns = {
        "bams": "bam",
        "bedMethyl": "bed",
        "openChromatin": "bed",
        "annot": ["bam", "csv", "tsv", "gtf"]
    }
    print("Scanning for files to inventory...")
    for subfolder, extensions in search_patterns.items():
        target_dir = os.path.join(directory, subfolder)
        if not os.path.exists(target_dir):
            continue
        if not isinstance(extensions, list):
            extensions = [extensions]
        for ext in extensions:
            pattern = os.path.join(target_dir, f'*.{ext}')
            for filepath_str in glob.glob(pattern):
                filepath = Path(filepath_str)
                print(f"  - Processing inventory for: {filepath.name}")
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

def write_inventory_report(data, output_file):
    """Writes the collected file inventory data to a TSV file."""
    with open(output_file, 'w', newline='') as csvfile:
        fieldnames = ['filename', 'extension', 'size_bytes', 'path', 'mod_time', 'md5sum']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerows(data)

# ==============================================================================
# SECTION 2: QC STATISTICS FUNCTIONS
# ==============================================================================

def read_fastq_lengths_and_q(path):
    """Processes a FASTQ file to get read lengths and quality scores."""
    lengths, n_reads, total_bases, sum_q = [], 0, 0, 0.0
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="ignore") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1: # Sequence line
                seq = line.strip()
                L = len(seq)
                lengths.append(L)
                total_bases += L
            elif i % 4 == 3: # Quality line
                qual = line.strip()
                s = sum((ord(c) - 33) for c in qual)
                sum_q += (s / len(qual)) if len(qual) > 0 else 0.0
                n_reads += 1
    return lengths, n_reads, total_bases, (sum_q / n_reads if n_reads else None)

def read_bam_lengths(path):
    """Processes a BAM file to get read lengths and mapping stats."""
    if pysam is None:
        print("[ERROR] 'pysam' library is not installed. Cannot process BAM files for QC.", file=sys.stderr)
        return [], 0, 0, 0
    lengths, n_reads, total_bases, mapped = [], 0, 0, 0
    with pysam.AlignmentFile(path, "rb", check_sq=False) as bam:
        for r in bam.fetch(until_eof=True):
            if r.is_secondary or r.is_supplementary:
                continue
            L = r.query_length or 0
            lengths.append(L)
            n_reads += 1
            total_bases += L
            if not r.is_unmapped:
                mapped += 1
    return lengths, n_reads, total_bases, mapped

def percentile(sorted_list, q):
    """Calculates the q-th percentile of a sorted list."""
    if not sorted_list: return None
    k = (len(sorted_list) - 1) * (q / 100.0)
    f, c = int(k), min(int(k) + 1, len(sorted_list) - 1)
    if f == c: return float(sorted_list[f])
    return sorted_list[f] * (c - k) + sorted_list[c] * (k - f)

def n50(lengths):
    """Calculates the N50 for a list of read lengths."""
    tot = sum(lengths)
    half = tot / 2
    for L in sorted(lengths, reverse=True):
        half -= L
        if half <= 0: return L
    return None

def compute_qc_stats(lengths):
    """Computes a dictionary of statistics from a list of lengths."""
    if not lengths: return {}
    s = sorted(lengths); n = len(s); tot = sum(s)
    return {
        "mean_len": tot / n, "median_len": median(s), "N50_len": n50(s),
        "p10_len": percentile(s, 10), "p25_len": percentile(s, 25),
        "p75_len": percentile(s, 75), "p90_len": percentile(s, 90),
        "min_len": s[0], "max_len": s[-1],
        "ge_1kb_pct": 100 * sum(L >= 1000 for L in s) / n,
        "ge_5kb_pct": 100 * sum(L >= 5000 for L in s) / n,
        "ge_10kb_pct": 100 * sum(L >= 10000 for L in s) / n
    }

def write_qc_row(outf, sample, file_type, filename, n_reads, total_bases, meanQ, mapped, mapped_pct, stats):
    """Writes a single row to the QC statistics CSV file."""
    yieldGb = total_bases / 1e9 if total_bases else 0
    vals = [
        sample, file_type, filename, n_reads, total_bases, f"{yieldGb:.3f}",
        f"{stats.get('mean_len', 0):.2f}", f"{stats.get('median_len', 0):.0f}", f"{stats.get('N50_len', 0):.0f}",
        f"{stats.get('p10_len', 0):.0f}", f"{stats.get('p25_len', 0):.0f}", f"{stats.get('p75_len', 0):.0f}", f"{stats.get('p90_len', 0):.0f}",
        f"{stats.get('min_len', 0):.0f}", f"{stats.get('max_len', 0):.0f}",
        f"{stats.get('ge_1kb_pct', 0):.2f}", f"{stats.get('ge_5kb_pct', 0):.2f}", f"{stats.get('ge_10kb_pct', 0):.2f}",
        f"{meanQ:.2f}" if meanQ is not None else "", mapped or "", f"{mapped_pct:.2f}" if mapped_pct is not None else ""
    ]
    outf.write(",".join(str(x) for x in vals) + "\n")


# ==============================================================================
# SECTION 3: MAIN EXECUTION LOGIC
# ==============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Generate a file inventory report and a QC statistics report for a Dogme run.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    # This script now requires arguments for BOTH reports.
    # The Nextflow script MUST be updated to provide these arguments.
    parser.add_argument("-i", "--input_dir", required=True, help="Main directory to scan for analysis files.")
    parser.add_argument("-o_inv", "--output_inventory", required=True, help="Output path for the file inventory report (.tsv).")
    parser.add_argument("-o_qc", "--output_qc", required=True, help="Output path for the QC statistics report (.csv).")
    parser.add_argument("-s", "--sample", default=None, help="Sample name for the QC report (inferred from filenames if not provided).")
    args = parser.parse_args()

    # --- Part 1: Generate File Inventory Report ---
    print("--- Generating File Inventory Report (Part 1/2) ---")
    inventory_data = collect_file_info(args.input_dir)
    if inventory_data:
        write_inventory_report(inventory_data, args.output_inventory)
        print(f"\n[SUCCESS] File inventory report written to: {args.output_inventory}")
    else:
        print("\n[WARNING] No files found for the inventory report.")

    # --- Part 2: Generate QC Statistics Report ---
    print("\n--- Generating QC Statistics Report (Part 2/2) ---")
    qc_header = [
        "sample", "file_type", "filename", "n_reads", "total_bases", "yield_Gb",
        "mean_len", "median_len", "N50_len", "p10_len", "p25_len", "p75_len", "p90_len",
        "min_len", "max_len", "reads_ge_1kb_pct", "reads_ge_5kb_pct", "reads_ge_10kb_pct",
        "mean_per_read_Q", "mapped_reads", "mapped_pct"
    ]
    with open(args.output_qc, "w", newline='') as outf:
        outf.write(",".join(qc_header) + "\n")
        # Process FASTQ files
        fqdir = os.path.join(args.input_dir, "fastqs")
        if os.path.isdir(fqdir):
            for f in os.listdir(fqdir):
                if f.endswith((".fastq", ".fastq.gz")):
                    print(f"  - Processing QC for FASTQ: {f}")
                    path = os.path.join(fqdir, f)
                    sample = args.sample or Path(f).name.replace(".fastq.gz", "").replace(".fastq", "")
                    lengths, n_reads, tot, meanQ = read_fastq_lengths_and_q(path)
                    stats = compute_qc_stats(lengths)
                    write_qc_row(outf, sample, "FASTQ", f, n_reads, tot, meanQ, None, None, stats)
        # Process BAM files
        bamdir = os.path.join(args.input_dir, "bams")
        if os.path.isdir(bamdir):
            for f in os.listdir(bamdir):
                if f.endswith(".bam") and "plus" not in f and "minus" not in f:
                    print(f"  - Processing QC for BAM: {f}")
                    path = os.path.join(bamdir, f)
                    sample = args.sample or Path(f).name.replace(".bam", "")
                    lengths, n_reads, tot, mapped = read_bam_lengths(path)
                    stats = compute_qc_stats(lengths)
                    mapped_pct = (100 * mapped / n_reads) if n_reads else None
                    ftype = "Unmapped BAM" if "unmapped" in f else "BAM"
                    write_qc_row(outf, sample, ftype, f, n_reads, tot, None, mapped, mapped_pct, stats)
        # Process Annotated BAMs
        andir = os.path.join(args.input_dir, "annot")
        if os.path.isdir(andir):
            for f in os.listdir(andir):
                if f.endswith(".annotated.bam"):
                    print(f"  - Processing QC for Annotated BAM: {f}")
                    path = os.path.join(andir, f)
                    sample = args.sample or Path(f).name.replace(".annotated.bam", "")
                    lengths, n_reads, tot, mapped = read_bam_lengths(path)
                    stats = compute_qc_stats(lengths)
                    mapped_pct = (100 * mapped / n_reads) if n_reads else None
                    write_qc_row(outf, sample, "Annotated BAM", f, n_reads, tot, None, mapped, mapped_pct, stats)

    print(f"\n[SUCCESS] QC statistics report written to: {args.output_qc}")
    print("\nAll reports generated successfully. ✨")

if __name__ == "__main__":
    main()