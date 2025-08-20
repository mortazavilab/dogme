#!/usr/bin/env python3
import argparse
import os
import sys
import re
from typing import Tuple, Dict, Set, List
import pysam

# ==============================================================================
# Novel / Known classification (unchanged)
# ==============================================================================

DEFAULT_NOVEL_TT = {"nnc", "nic", "ism", "ism-j", "ism-n", "novel", "fusion"}
KNOWN_TT_BLOCK = {"known", "antisense"}

KNOWN_TX_RE = re.compile(r"^(ENST|ENSG|NM_|NR_)", re.IGNORECASE)
NOVEL_TX_RE_STRICT = re.compile(r"^NOVEL", re.IGNORECASE)
ALREADY_UNIFIED_TXRE = re.compile(r"^CONST_|^TALONT_", re.IGNORECASE)
ALREADY_UNIFIED_GXRE = re.compile(r"^CONSG_|^TALONG_", re.IGNORECASE)

def parse_novel_tt(arg: str) -> Set[str]:
    if not arg:
        return set(DEFAULT_NOVEL_TT)
    return {x.strip().lower() for x in arg.split(",") if x.strip()}

def is_novel_read(r: pysam.AlignedSegment, novel_tt_set: Set[str]) -> bool:
    if r.has_tag("TX") and ALREADY_UNIFIED_TXRE.match(str(r.get_tag("TX"))):
        return False
    if r.has_tag("GX") and ALREADY_UNIFIED_GXRE.match(str(r.get_tag("GX"))):
        return False
    if r.has_tag("TX") and KNOWN_TX_RE.match(str(r.get_tag("TX"))):
        return False
    if r.has_tag("TT"):
        tt = str(r.get_tag("TT")).strip().lower()
        if tt in KNOWN_TT_BLOCK:
            return False
        return tt in novel_tt_set
    if r.has_tag("TX") and NOVEL_TX_RE_STRICT.match(str(r.get_tag("TX"))):
        return True
    if r.has_tag("GX") and str(r.get_tag("GX")).strip().lower().startswith("novel"):
        return True
    return False

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

def get_exon_blocks(read: pysam.AlignedSegment) -> List[Tuple[int, int]]:
    return read.get_blocks()

# ==============================================================================
# MODIFIED SECTION
# ==============================================================================

def get_unique_sample_name(path: str, all_paths: List[str]) -> str:
    """Creates a unique, clean sample name from a file path."""
    name = os.path.basename(path).replace(".bam", "")
    # If basenames are not unique, prepend the parent directory name
    basenames = [os.path.basename(p) for p in all_paths]
    if basenames.count(os.path.basename(path)) > 1:
        parent_dir = os.path.basename(os.path.dirname(path))
        name = f"{parent_dir}_{name}"
    return name

def collect_and_assign_ids(
    bam_files: List[str],
    sample_names: Dict[str, str], # Pass in the unique names
    id_tag: str,
    gene_tag: str,
    gene_prefix: str,
    tx_prefix: str,
    novel_tt_set: Set[str]
) -> Tuple[Dict, Dict, Dict]:
    master_structures = {}
    old_tx_to_struct = {}
    gene_counter, transcript_counter = 1, 1

    for bam_path in bam_files:
        unique_name = sample_names[bam_path] # Use the pre-generated unique name
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for r in bam:
                if r.is_unmapped or not r.cigartuples: continue
                key = splice_key(r)
                if r.has_tag(id_tag):
                    old_tx_to_struct[str(r.get_tag(id_tag))] = key
                if key not in master_structures:
                    master_structures[key] = {
                        "gene_id": None, "transcript_id": None,
                        "is_novel": True,
                        "exon_blocks": get_exon_blocks(r),
                        # Initialize counts for all unique sample names
                        "read_counts": {s_name: 0 for s_name in sample_names.values()}
                    }
                # Update read count for the current UNIQUE sample
                master_structures[key]["read_counts"][unique_name] += 1
                if not is_novel_read(r, novel_tt_set):
                    if master_structures[key]["is_novel"]:
                        master_structures[key]["is_novel"] = False
                        master_structures[key]["gene_id"] = str(r.get_tag(gene_tag)) if r.has_tag(gene_tag) else "UnknownGene"
                        master_structures[key]["transcript_id"] = str(r.get_tag(id_tag)) if r.has_tag(id_tag) else "UnknownTranscript"

    for struct, data in master_structures.items():
        if data["is_novel"]:
            data["gene_id"] = f"{gene_prefix}{gene_counter}"
            data["transcript_id"] = f"{tx_prefix}{transcript_counter}"
            gene_counter += 1
            transcript_counter += 1

    struct_to_new_id = {k: (v["gene_id"], v["transcript_id"]) for k, v in master_structures.items()}
    return master_structures, struct_to_new_id, old_tx_to_struct

# ==============================================================================
# Unchanged sections...
# ==============================================================================

def rewrite_bam(in_bam: str, out_bam: str, tx_map: Dict, struct_map: Dict, id_tag: str, gene_tag: str, novel_tt_set: Set[str]) -> None:
    changed_count = 0; total_count = 0
    with pysam.AlignmentFile(in_bam, "rb") as ib, pysam.AlignmentFile(out_bam, "wb", header=ib.header) as ob:
        for r in ib:
            total_count += 1
            if not is_novel_read(r, novel_tt_set):
                ob.write(r)
                continue
            new_gx, new_tx = None, None
            if r.has_tag(id_tag):
                tx_id = str(r.get_tag(id_tag))
                if tx_id in tx_map: new_gx, new_tx = tx_map[tx_id]
            if new_gx is None:
                struct = splice_key(r)
                if struct in struct_map: new_gx, new_tx = struct_map[struct]
            if new_gx and new_tx:
                r.set_tag(gene_tag, new_gx, value_type='Z')
                r.set_tag(id_tag, new_tx, value_type='Z')
                changed_count += 1
            ob.write(r)
    try: pysam.index(out_bam)
    except Exception as e: print(f"[WARN] Could not index {out_bam}: {e}", file=sys.stderr)
    print(f"  - Rewrote {os.path.basename(in_bam)}: changed IDs for {changed_count} / {total_count} reads.")

def main():
    ap = argparse.ArgumentParser(description="Generate TALON-compatible files and reconciled BAMs, preserving known annotations.")
    ap.add_argument("--bams", required=True, nargs='+', help="List of input BAM files.")
    ap.add_argument("--out_prefix", required=True, help="Prefix for GTF and abundance files.")
    ap.add_argument("--outdir", required=True, help="Directory to save all output files.")
    ap.add_argument("--gene_prefix", default="TALONG_", help="Consolidated novel gene ID prefix.")
    ap.add_argument("--tx_prefix", default="TALONT_", help="Consolidated novel transcript ID prefix.")
    ap.add_argument("--id_tag", default="TX", help="BAM tag for transcript ID.")
    ap.add_argument("--gene_tag", default="GX", help="BAM tag for gene ID.")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    novel_tt_set = parse_novel_tt("")

    # === MODIFIED SECTION: Create unique names upfront ===
    sample_names_map = {path: get_unique_sample_name(path, args.bams) for path in args.bams}
    unique_sample_names = list(sample_names_map.values())

    print("=== Step 1 & 2: Collecting structures and assigning IDs... ===")
    all_data, struct_to_new_id, old_tx_to_struct = collect_and_assign_ids(
        args.bams, sample_names_map, args.id_tag, args.gene_tag, args.gene_prefix, args.tx_prefix, novel_tt_set
    )
    print(f"Processed {len(all_data)} unique transcript structures across {len(args.bams)} samples.")

    old_tx_to_new_id = {tx: struct_to_new_id.get(struct) for tx, struct in old_tx_to_struct.items() if struct_to_new_id.get(struct)}
    print(f"Mapped {len(old_tx_to_new_id)} original transcript IDs to final IDs.")

    print("\n=== Step 3: Rewriting BAM files (novel reads only)... ===")
    for bam_file in args.bams:
        base_name = os.path.basename(bam_file)
        out_bam_path = os.path.join(args.outdir, f"{os.path.splitext(base_name)[0]}.reconciled.bam")
        rewrite_bam(bam_file, out_bam_path, old_tx_to_new_id, struct_to_new_id, args.id_tag, args.gene_tag, novel_tt_set)

    print("\n=== Step 4: Generating TALON files... ===")
    gtf_file = os.path.join(args.outdir, f"{args.out_prefix}.gtf")
    with open(gtf_file, 'w') as f:
        for struct, data in all_data.items():
            gene_id, tx_id = data["gene_id"], data["transcript_id"]
            chrom, strand, _ = struct
            min_start = min(block[0] for block in data["exon_blocks"]) + 1
            max_end = max(block[1] for block in data["exon_blocks"])
            attributes = f'gene_id "{gene_id}"; transcript_id "{tx_id}";'
            f.write(f"{chrom}\tTALON\ttranscript\t{min_start}\t{max_end}\t.\t{strand}\t.\t{attributes}\n")
            for exon_start, exon_end in data["exon_blocks"]:
                f.write(f"{chrom}\tTALON\texon\t{exon_start + 1}\t{exon_end}\t.\t{strand}\t.\t{attributes}\n")

    abundance_file = os.path.join(args.outdir, f"{args.out_prefix}_abundance.tsv")
    with open(abundance_file, 'w') as f:
        # Use the unique_sample_names list for the header
        header = ["gene_ID", "transcript_ID", "annot_gene_id", "annot_transcript_id",
                  "annot_gene_name", "annot_transcript_name", "n_exons",
                  "transcript_length", "gene_novelty", "transcript_novelty",
                  "ISM_subtype"] + unique_sample_names
        f.write('\t'.join(header) + '\n')

        for struct, data in all_data.items():
            gene_id, tx_id = data["gene_id"], data["transcript_id"]
            novelty_str = "Novel" if data["is_novel"] else "Known"
            annot_gene_id = "NA" if data["is_novel"] else gene_id
            annot_tx_id = "NA" if data["is_novel"] else tx_id
            row = [gene_id, tx_id, annot_gene_id, annot_tx_id, "NA", "NA", str(len(data['exon_blocks'])), "NA", novelty_str, novelty_str, "NA"]
            # Loop through the unique names to get the correct counts
            for s_name in unique_sample_names:
                row.append(str(data["read_counts"].get(s_name, 0)))
            f.write('\t'.join(row) + '\n')

    print("\n=== ✅ All tasks complete! ===")
    print(f"Outputs are in: {args.outdir}")

if __name__ == "__main__":
    main()