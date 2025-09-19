#!/usr/bin/env python3
import argparse
import os
import sys
import re
from typing import Tuple, Dict, Set, List
from collections import Counter
import pysam
import concurrent.futures

__version__ = "1.0.4"

# ==============================================================================
# GTF Parsing
# ==============================================================================
def parse_gtf_for_ids_and_names(gtf_file: str) -> Tuple[Set[str], Set[str], Dict[str, str], Dict[str, str], Dict[str, List[str]], Dict[str, List[Tuple[int, int]]]]:
    if not gtf_file or not os.path.exists(gtf_file):
        print("[WARN] No valid annotation GTF provided...", file=sys.stderr)
        return set(), set(), {}, {}, {}, {}

    print(f"=== Parsing reference annotation for known IDs and names: {os.path.basename(gtf_file)}... ===")
    known_gene_ids, known_transcript_ids = set(), set()
    gene_id_to_name, transcript_id_to_name = {}, {}
    transcript_gtf_lines = {}
    transcript_exon_blocks = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'): continue
            parts = line.strip().split('\t')
            if len(parts) < 9: continue
            attributes_str = parts[8]
            gene_id_match = re.search(r'gene_id\s+"([^"]+)"', attributes_str)
            tx_id_match = re.search(r'transcript_id\s+"([^"]+)"', attributes_str)
            gene_name_match = re.search(r'gene_name\s+"([^"]+)"', attributes_str)
            tx_name_match = re.search(r'transcript_name\s+"([^"]+)"', attributes_str)
            if gene_id_match:
                base_gene_id = gene_id_match.group(1).split('.')[0]
                known_gene_ids.add(base_gene_id)
                if gene_name_match and base_gene_id not in gene_id_to_name:
                    gene_id_to_name[base_gene_id] = gene_name_match.group(1)
            if tx_id_match:
                base_tx_id = tx_id_match.group(1).split('.')[0]
                known_transcript_ids.add(base_tx_id)
                if tx_name_match and base_tx_id not in transcript_id_to_name:
                    transcript_id_to_name[base_tx_id] = tx_name_match.group(1)
                # Store all lines for this transcript
                if base_tx_id not in transcript_gtf_lines:
                    transcript_gtf_lines[base_tx_id] = []
                transcript_gtf_lines[base_tx_id].append(line.rstrip('\n'))
                # Collect exon blocks for this transcript
                if parts[2] == "exon":
                    start = int(parts[3]) - 1
                    end = int(parts[4])
                    if base_tx_id not in transcript_exon_blocks:
                        transcript_exon_blocks[base_tx_id] = []
                    transcript_exon_blocks[base_tx_id].append((start, end))
    print(f"Found {len(known_gene_ids)} known gene IDs and {len(known_transcript_ids)} known transcript IDs.")
    print(f"Found names for {len(gene_id_to_name)} genes and {len(transcript_id_to_name)} transcripts.")
    return known_gene_ids, known_transcript_ids, gene_id_to_name, transcript_id_to_name, transcript_gtf_lines, transcript_exon_blocks

# ==============================================================================
# Splicing and Exon Definitions
# ==============================================================================
def splice_key(read: pysam.AlignedSegment) -> Tuple[str, str, Tuple[str, ...]]:
    chrom = read.reference_name
    strand = '-' if read.is_reverse else '+'
    juncs = []
    ref_pos = read.reference_start
    ref_consume = {pysam.CMATCH, pysam.CDEL, pysam.CREF_SKIP, pysam.CEQUAL, pysam.CDIFF}
    if read.cigartuples:
        for op, length in read.cigartuples:
            if op == pysam.CREF_SKIP:
                juncs.append(f"{chrom}:{ref_pos+1}-{ref_pos+length}")
            if op in ref_consume:
                ref_pos += length
    return (chrom, strand, tuple(sorted(juncs)))

def get_exon_blocks(read: pysam.AlignedSegment, merge_distance: int = 5) -> List[Tuple[int, int]]:
    """
    Returns a list of exon blocks for the read, merging blocks that are <= merge_distance apart.
    """
    blocks = read.get_blocks()
    if not blocks:
        return []
    merged = [list(blocks[0])]
    for start, end in blocks[1:]:
        prev_end = merged[-1][1]
        if start - prev_end <= merge_distance:
            merged[-1][1] = end
        else:
            merged.append([start, end])
    return [tuple(b) for b in merged]

# ==============================================================================
# Multithreading Logic
# ==============================================================================
def get_unique_sample_name(path: str, all_paths: List[str]) -> str:
    name = os.path.basename(path).replace(".bam", "")
    basenames = [os.path.basename(p) for p in all_paths]
    if basenames.count(os.path.basename(path)) > 1:
        parent_dir = os.path.basename(os.path.dirname(path))
        name = f"{parent_dir}_{name}"
    return name

def _collect_reads_worker(bam_path: str, id_tag: str, gene_tag: str, merge_distance: int) -> Dict:
    local_reads = {}
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for r in bam:
            if r.is_unmapped or not r.cigartuples or len(r.cigartuples) < 3: continue
            key = splice_key(r)
            if not key[2]: continue
            if key not in local_reads:
                local_reads[key] = {
                    "lengths": [],
                    "rep_read": {
                        "tx_id": str(r.get_tag(id_tag)) if r.has_tag(id_tag) else None,
                        "gx_id": str(r.get_tag(gene_tag)) if r.has_tag(gene_tag) else None,
                        "tt_tag": str(r.get_tag('TT')) if r.has_tag('TT') else "Novel",
                        "exon_blocks": get_exon_blocks(r, merge_distance)
                    }
                }
            local_reads[key]["lengths"].append(r.query_alignment_length)
    return (bam_path, local_reads)

def collect_and_assign_ids(
    bam_files: List[str], sample_names: Dict[str, str], known_gene_ids: Set[str],
    known_tx_ids: Set[str], id_tag: str, gene_tag: str, gene_prefix: str,
    tx_prefix: str, threads: int, merge_distance: int
) -> Dict:
    all_unique_samples = list(sample_names.values())
    master_structures = {}

    with concurrent.futures.ThreadPoolExecutor(max_workers=threads) as executor:
        future_to_bam = {
            executor.submit(_collect_reads_worker, bam_path, id_tag, gene_tag, merge_distance): bam_path
            for bam_path in bam_files
        }
        for future in concurrent.futures.as_completed(future_to_bam):
            try:
                bam_path, local_reads = future.result()
                unique_name = sample_names[bam_path]
                for key, data in local_reads.items():
                    if key not in master_structures:
                        master_structures[key] = {"read_counts": {s: 0 for s in all_unique_samples}, "rep_read": data["rep_read"], "lengths": []}
                    master_structures[key]["read_counts"][unique_name] += len(data["lengths"])
                    master_structures[key]["lengths"].extend(data["lengths"])
            except Exception as exc:
                print(f"{future_to_bam[future]} generated an exception: {exc}", file=sys.stderr)

    gene_counter, transcript_counter = 1, 1
    for data in master_structures.values():
        rep = data["rep_read"]
        if data["lengths"]: data["transcript_length"] = round(sum(data["lengths"]) / len(data["lengths"]))
        else: data["transcript_length"] = 0
        data["exon_blocks"] = rep["exon_blocks"]
        tt_tag = rep.get("tt_tag", "Novel").upper()
        data["tt_tag"] = tt_tag
        read_tx_id = rep["tx_id"]
        read_gene_id = rep["gx_id"]

        if tt_tag in ("KNOWN", "ISM"):
            data["transcript_id"] = read_tx_id
            data["gene_id"] = read_gene_id
            data["is_novel_gene"] = False
            data["is_novel_transcript"] = (tt_tag == "ISM")
        else:
            data["is_novel_transcript"] = True
            if read_gene_id and read_gene_id in known_gene_ids:
                data["gene_id"] = read_gene_id
                data["is_novel_gene"] = False
            else:
                data["gene_id"] = None
                data["is_novel_gene"] = True
            if sum(data["read_counts"].values()) == 1:
                data["transcript_id"] = "solo"
            else:
                data["transcript_id"] = None

    for data in master_structures.values():
        if data["transcript_id"] is None:
            data["transcript_id"] = f"{tx_prefix}{transcript_counter}"
            transcript_counter += 1
        if data["gene_id"] is None:
            data["gene_id"] = f"{gene_prefix}{gene_counter}"
            gene_counter += 1
    return master_structures

# ==============================================================================
# ISM Consolidation Logic
# ==============================================================================
def consolidate_transcript_variants(raw_structures: Dict, transcript_exon_blocks: Dict[str, List[Tuple[int, int]]]) -> Tuple[Dict, Dict]:
    """
    Consolidates multiple variants of the same transcript ID based on TT tag.
    All 'ISM' variants for a given ID are merged into one.
    All 'Known' variants for a given ID are merged into one, and the representative's exon_blocks are set to match the annotation GTF.
    """
    final_structures = {}
    key_remap = {}
    representatives = {}

    for splice_key, data in raw_structures.items():
        tt_tag = data.get("tt_tag")
        tx_id = data.get("transcript_id")

        if tt_tag == "KNOWN":
            group_key = (tx_id, tt_tag)
            if group_key not in representatives:
                representatives[group_key] = splice_key
                # Overwrite exon_blocks and transcript_length to match annotation GTF
                if tx_id in transcript_exon_blocks:
                    data["exon_blocks"] = sorted(transcript_exon_blocks[tx_id], key=lambda x: x[0])
                    if data["exon_blocks"]:
                        data["transcript_length"] = sum(e - s for s, e in data["exon_blocks"])
                final_structures[splice_key] = data
            else:
                rep_key = representatives[group_key]
                for sample, count in data["read_counts"].items():
                    final_structures[rep_key]["read_counts"][sample] += count
                final_structures[rep_key]["lengths"].extend(data["lengths"])
                key_remap[splice_key] = rep_key
        elif tt_tag == "ISM":
            group_key = (tx_id, tt_tag)
            if group_key not in representatives:
                representatives[group_key] = splice_key
                final_structures[splice_key] = data
            else:
                rep_key = representatives[group_key]
                for sample, count in data["read_counts"].items():
                    final_structures[rep_key]["read_counts"][sample] += count
                final_structures[rep_key]["lengths"].extend(data["lengths"])
                key_remap[splice_key] = rep_key
        else:
            final_structures[splice_key] = data

    # Recalculate mean length for the merged groups
    for rep_key in representatives.values():
        if rep_key in final_structures:
            lengths = final_structures[rep_key]["lengths"]
            if lengths:
                final_structures[rep_key]["transcript_length"] = round(sum(lengths) / len(lengths))
            # For KNOWN, ensure transcript_length matches annotation if available
            data = final_structures[rep_key]
            if data.get("tt_tag") == "KNOWN":
                tx_id = data.get("transcript_id")
                if tx_id in transcript_exon_blocks:
                    data["exon_blocks"] = sorted(transcript_exon_blocks[tx_id], key=lambda x: x[0])
                    if data["exon_blocks"]:
                        data["transcript_length"] = sum(e - s for s, e in data["exon_blocks"])

    print(f"Consolidated {len(key_remap)} redundant KNOWN/ISM structures into representative entries.")
    return final_structures, key_remap

# ==============================================================================
# Rewriting and Output Generation
# ==============================================================================
def rewrite_bam(in_bam: str, out_bam: str, struct_to_new_id: Dict, id_tag: str, gene_tag: str) -> None:
    changed_count, total_count = 0, 0
    mappings = set()
    # Read input BAM and header
    with pysam.AlignmentFile(in_bam, "rb") as ib:
        header = ib.header.to_dict()
        # --- Add/merge PG (program) and CO (comment) lines for provenance ---
        # Copy existing PG/CO lines if present
        pg_lines = header.get('PG', [])
        co_lines = header.get('CO', [])
        # Add a new PG line for this program
        new_pg = {
            'ID': 'reconcileBams',
            'PN': 'reconcileBams.py',
            'VN': __version__,
            'CL': ' '.join(sys.argv)
        }
        pg_lines.append(new_pg)
        header['PG'] = pg_lines
        # Add a comment line for this run
        co_lines.append(f"reconcileBams.py version {__version__} command: {' '.join(sys.argv)}")
        header['CO'] = co_lines

        with pysam.AlignmentFile(out_bam, "wb", header=header) as ob:
            for r in ib:
                total_count += 1
                key = splice_key(r)
                if key in struct_to_new_id:
                    new_gx, new_tx = struct_to_new_id[key]
                    current_tx = str(r.get_tag(id_tag)) if r.has_tag(id_tag) else "NA"
                    current_gx = str(r.get_tag(gene_tag)) if r.has_tag(gene_tag) else "NA"
                    if new_tx != current_tx or new_gx != current_gx:
                        if new_gx: r.set_tag(gene_tag, new_gx, value_type='Z')
                        if new_tx: r.set_tag(id_tag, new_tx, value_type='Z')
                        changed_count += 1
                        if new_tx != current_tx:
                            mappings.add((current_tx, new_tx))
                ob.write(r)
    mapping_file_path = out_bam.replace(".reconciled.bam", ".mapping.tsv")
    with open(mapping_file_path, "w") as f_map:
        f_map.write("original_tx_id\tfinal_tx_id\n")
        for orig_id, final_id in sorted(list(mappings)):
            f_map.write(f"{orig_id}\t{final_id}\n")
    try: pysam.index(out_bam)
    except Exception as e: print(f"[WARN] Could not index {out_bam}: {e}", file=sys.stderr)
    print(f"  - Finished rewriting {os.path.basename(in_bam)}: changed IDs for {changed_count} / {total_count} reads. Mapping file created.")

def main():
    ap = argparse.ArgumentParser(description="Generate reconciled BAM and annotation files.")
    ap.add_argument("--bams", required=True, nargs='+', help="List of input BAM files.")
    ap.add_argument("--annotation", required=True, help="Reference annotation GTF file.")
    ap.add_argument("--out_prefix", required=True, help="Prefix for GTF and abundance files.")
    ap.add_argument("--outdir", required=True, help="Directory to save all output files.")
    ap.add_argument("--gene_prefix", default="CONSG", help="Consolidated novel gene ID prefix.")
    ap.add_argument("--tx_prefix", default="CONST", help="Consolidated novel transcript ID prefix.")
    ap.add_argument("--id_tag", default="TX", help="BAM tag for transcript ID.")
    ap.add_argument("--gene_tag", default="GX", help="BAM tag for gene ID.")
    ap.add_argument("--threads", type=int, default=os.cpu_count(), help="Number of threads to use.")
    ap.add_argument("--exon_merge_distance", type=int, default=5, help="Merge exon blocks that are <= this many nucleotides apart (default: 5).")
    args = ap.parse_args()

    print(f"reconcileBams.py version {__version__}")

    os.makedirs(args.outdir, exist_ok=True)
    known_gene_ids, known_tx_ids, gene_id_to_name, transcript_id_to_name, transcript_gtf_lines, transcript_exon_blocks = parse_gtf_for_ids_and_names(args.annotation)
    sample_names_map = {path: get_unique_sample_name(path, args.bams) for path in args.bams}

    print(f"\n=== Step 1: Collecting structures and assigning IDs (using {args.threads} threads)... ===")
    raw_data = collect_and_assign_ids(
        args.bams, sample_names_map, known_gene_ids, known_tx_ids,
        args.id_tag, args.gene_tag, args.gene_prefix, args.tx_prefix, args.threads,
        args.exon_merge_distance
    )
    print(f"Initially processed {len(raw_data)} unique transcript structures.")

    print("\n=== Step 1.5: Consolidating KNOWN/ISM variants... ===")
    final_data, key_remap = consolidate_transcript_variants(raw_data, transcript_exon_blocks)
    print(f"Final dataset contains {len(final_data)} unique transcripts after consolidation.")

    # --- Summary of Findings ---
    print("\n=== Summary of Findings ===")
    summary_lines = []
    category_counts = Counter(data.get('tt_tag', 'UNKNOWN') for data in final_data.values() if data.get('transcript_id') != 'solo')
    solo_count = sum(1 for data in final_data.values() if data.get('transcript_id') == 'solo')

    # Only count truly novel transcript IDs (not ISM/KNOWN) for summary
    novel_transcript_ids = set()
    for data in final_data.values():
        tt = data.get('tt_tag', '').upper()
        tid = data.get('transcript_id')
        if data.get('is_novel_transcript', False) and tid and tid != 'solo' and tt not in ("KNOWN", "ISM"):
            novel_transcript_ids.add(tid)
    # Count novel gene IDs as before
    novel_gene_ids = set()
    for data in final_data.values():
        if data.get('is_novel_gene', False) and data.get('gene_id'):
            novel_gene_ids.add(data['gene_id'])

    for category, count in sorted(category_counts.items()):
        line = f"  - {category.capitalize()} transcripts:{' ':<27} {count}"
        print(line)
        summary_lines.append(line)
    if solo_count > 0:
        line = f"  - {'Single-read solo transcripts (will be filtered):':<45} {solo_count}"
        print(line)
        summary_lines.append(line)
    line = f"  - Number of genes with {args.gene_prefix} IDs: {len(novel_gene_ids)}"
    print(line)
    summary_lines.append(line)
    line = f"  - Number of transcripts with {args.tx_prefix} IDs: {len(novel_transcript_ids)}"
    print(line)
    summary_lines.append(line)

    # Write summary to report file
    report_file = os.path.join(args.outdir, f"{args.out_prefix}_summary.txt")
    with open(report_file, "w") as f:
        f.write(f"reconcileBams.py version {__version__}\n")
        f.write("=== Summary of Findings ===\n")
        for line in summary_lines:
            f.write(line + "\n")

    struct_to_new_id = {k: (v["gene_id"], v["transcript_id"]) for k, v in final_data.items()}
    for merged_key, representative_key in key_remap.items():
        if representative_key in struct_to_new_id:
            struct_to_new_id[merged_key] = struct_to_new_id[representative_key]

    print(f"\n=== Step 2: Rewriting BAM files (using {args.threads} threads)... ===")
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(
            rewrite_bam, bam_file, os.path.join(args.outdir, f"{os.path.splitext(os.path.basename(bam_file))[0]}.reconciled.bam"),
            struct_to_new_id, args.id_tag, args.gene_tag
        ) for bam_file in args.bams]
        concurrent.futures.wait(futures)

    print("\n=== Step 3: Generating DOGME annotation files... ===")
    gtf_file = os.path.join(args.outdir, f"{args.out_prefix}.gtf")
    abundance_file = os.path.join(args.outdir, f"{args.out_prefix}_abundance.tsv")
    output_data = {k: v for k, v in final_data.items() if v["transcript_id"] != "solo"}
    print(f"Writing {len(output_data)} unique transcripts to output files.")

    with open(gtf_file, 'w') as f:
        for struct, data in output_data.items():
            gene_id, tx_id = data["gene_id"], data["transcript_id"]
            tt_tag = data.get("tt_tag", "").upper()
            if tt_tag == "KNOWN" and tx_id in transcript_gtf_lines:
                # MODIFICATION: Rewrite annotation lines to strip suffixes and add transcript_type
                for line in transcript_gtf_lines[tx_id]:
                    parts = line.strip().split('\t')
                    attributes_str = parts[8]
                    # Replace IDs with the suffix-less versions
                    attributes_str = re.sub(r'gene_id\s+"[^"]+"', f'gene_id "{gene_id}"', attributes_str)
                    attributes_str = re.sub(r'transcript_id\s+"[^"]+"', f'transcript_id "{tx_id}"', attributes_str)
                    # Add the new transcript_type attribute
                    attributes_str = f'{attributes_str.rstrip()} transcript_type "KNOWN";'
                    parts[8] = attributes_str
                    f.write('\t'.join(parts) + '\n')
                continue

            chrom, strand, _ = struct
            min_start = min(block[0] for block in data["exon_blocks"]) + 1
            max_end = max(block[1] for block in data["exon_blocks"])
            attributes = f'gene_id "{gene_id}"; transcript_id "{tx_id}"; transcript_type "{tt_tag}";'
            f.write(f"{chrom}\tDOGME\ttranscript\t{min_start}\t{max_end}\t.\t{strand}\t.\t{attributes}\n")
            for exon_start, exon_end in sorted(data["exon_blocks"], key=lambda x: x[0]):
                f.write(f"{chrom}\tDOGME\texon\t{exon_start + 1}\t{exon_end}\t.\t{strand}\t.\t{attributes}\n")

    with open(abundance_file, 'w') as f:
        unique_sample_names = list(sample_names_map.values())
        header = ["gene_ID", "transcript_ID", "annot_gene_id", "annot_transcript_id",
                  "annot_gene_name", "annot_transcript_name", "n_exons",
                  "transcript_length", "gene_novelty", "transcript_novelty",
                  "ISM_subtype"] + unique_sample_names
        f.write('\t'.join(header) + '\n')
        for struct, data in output_data.items():
            gene_id, tx_id = data["gene_id"], data["transcript_id"]
            annot_gene_id, annot_tx_id = gene_id, tx_id
            gene_novelty_str = "Known" if not data["is_novel_gene"] else "Novel"
            tx_novelty_str = data.get("tt_tag", "Novel")
            tx_len = str(data.get("transcript_length", "NA"))
            gene_name = gene_id_to_name.get(gene_id, "NA")
            tx_name = transcript_id_to_name.get(tx_id, "NA")
            row = [gene_id, tx_id, annot_gene_id, annot_tx_id, gene_name, tx_name, str(len(data['exon_blocks'])), tx_len, gene_novelty_str, tx_novelty_str, "NA"]
            for s_name in unique_sample_names:
                row.append(str(data["read_counts"].get(s_name, 0)))
            f.write('\t'.join(row) + '\n')

    print("\n=== ✅ All tasks complete! ===")
    print(f"Outputs are in: {args.outdir}")

if __name__ == "__main__":
    main()