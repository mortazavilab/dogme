#!/usr/bin/env python3

import pysam
import argparse
from intervaltree import Interval, IntervalTree
from collections import defaultdict
import sys
import os
import time

def parse_gtf_to_annotations(gtf_file_path: str):
    """
    Parses a GTF file to build data structures for genes and transcripts.

    Args:
        gtf_file_path: Path to the GTF annotation file.

    Returns:
        A tuple containing two dictionaries:
        - gene_annotations: Info about each gene, including its exons.
          {gene_id: {'chrom': str, 'strand': str, 'exons': set((start, end), ...)}}
        - transcript_annotations: Info about each transcript's exact exon structure.
          {transcript_id: {'gene_id': str, 'chrom': str, 'strand': str,
                           'exons': tuple(sorted((start, end), ...))}}
    """
    print(f"[{time.ctime()}] Parsing GTF file: {gtf_file_path}...")
    gene_annotations = defaultdict(lambda: {'exons': set(), 'chrom': None, 'strand': None})
    transcript_annotations = defaultdict(lambda: {'exons': set(), 'gene_id': None, 'chrom': None, 'strand': None})

    with open(gtf_file_path, 'r') as gtf:
        for line in gtf:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type != 'exon':
                continue

            chrom = fields[0]
            # GTF is 1-based, inclusive. Convert to 0-based, half-open for intervaltree/pysam.
            start = int(fields[3]) - 1
            end = int(fields[4])
            strand = fields[6]
            attributes = dict(item.strip().split(' ', 1) for item in fields[8].replace('"', '').split(';') if item)

            gene_id = attributes.get('gene_id')
            transcript_id = attributes.get('transcript_id')

            if not gene_id or not transcript_id:
                continue

            # Populate gene-level information
            gene_annotations[gene_id]['exons'].add((start, end))
            if not gene_annotations[gene_id]['chrom']:
                gene_annotations[gene_id]['chrom'] = chrom
                gene_annotations[gene_id]['strand'] = strand

            # Populate transcript-level information
            transcript_annotations[transcript_id]['exons'].add((start, end))
            if not transcript_annotations[transcript_id]['gene_id']:
                transcript_annotations[transcript_id]['gene_id'] = gene_id
                transcript_annotations[transcript_id]['chrom'] = chrom
                transcript_annotations[transcript_id]['strand'] = strand

    # Finalize transcript annotations by converting exon sets to sorted, hashable tuples
    final_transcript_annotations = {}
    for tx_id, data in transcript_annotations.items():
        final_transcript_annotations[tx_id] = {
            'gene_id': data['gene_id'],
            'chrom': data['chrom'],
            'strand': data['strand'],
            'exons': tuple(sorted(list(data['exons'])))
        }
        
    print(f"[{time.ctime()}] Parsed {len(gene_annotations)} genes and {len(final_transcript_annotations)} transcripts.")
    return dict(gene_annotations), final_transcript_annotations


def build_gene_interval_trees(gene_annotations: dict):
    """
    Builds interval trees for gene loci from parsed GTF data.

    Args:
        gene_annotations: The gene annotation dictionary from parse_gtf_to_annotations.

    Returns:
        A dictionary mapping chromosome names to IntervalTree objects.
        {chrom: IntervalTree([Interval(locus_start, locus_end, gene_id), ...])}
    """
    print(f"[{time.ctime()}] Building gene interval trees...")
    gene_interval_trees = defaultdict(IntervalTree)
    for gene_id, data in gene_annotations.items():
        if not data['exons']:
            continue
        
        # Define the gene locus as the min/max coordinates of all its exons
        locus_start = min(e[0] for e in data['exons'])
        locus_end = max(e[1] for e in data['exons'])
        chrom = data['chrom']
        
        gene_interval_trees[chrom].add(Interval(locus_start, locus_end, gene_id))
    
    print(f"[{time.ctime()}] Interval trees built for {len(gene_interval_trees)} chromosomes.")
    return gene_interval_trees

def get_read_strand(read: pysam.AlignedSegment) -> str:
    """Returns the strand of a read as '+' or '-'."""
    return '-' if read.is_reverse else '+'

def calculate_exonic_overlap(read_blocks: list, gene_exons: set) -> int:
    """
    Calculates the total length of overlap between a read's exons and a gene's exons.

    Args:
        read_blocks: A list of (start, end) tuples for the read's exons from get_blocks().
        gene_exons: A set of (start, end) tuples for the gene's exons.

    Returns:
        The total number of overlapping base pairs.
    """
    total_overlap = 0
    for r_start, r_end in read_blocks:
        for g_start, g_end in gene_exons:
            # Calculate overlap for each pair of exons
            overlap_start = max(r_start, g_start)
            overlap_end = min(r_end, g_end)
            overlap = max(0, overlap_end - overlap_start)
            total_overlap += overlap
    return total_overlap


def main():
    """Main execution function."""
    parser = argparse.ArgumentParser(
        description="Annotate a long-read RNA-seq BAM file with gene (gx) and transcript (tx) tags."
    )
    parser.add_argument("--bam", required=True, help="Input BAM file, sorted by coordinate.")
    parser.add_argument("--gtf", required=True, help="Input GTF annotation file.")
    parser.add_argument("--output", required=True, help="Output BAM file path.")
    args = parser.parse_args()

    # --- 1. Load and process annotations ---
    gene_annotations, transcript_annotations = parse_gtf_to_annotations(args.gtf)
    gene_interval_trees = build_gene_interval_trees(gene_annotations)
    
    # Create a reverse map from gene_id to a list of its transcript_ids
    gene_to_transcripts = defaultdict(list)
    for tx_id, data in transcript_annotations.items():
        gene_to_transcripts[data['gene_id']].append(tx_id)

    # --- 2. Setup for BAM processing ---
    print(f"[{time.ctime()}] Starting BAM file processing...")
    try:
        infile = pysam.AlignmentFile(args.bam, "rb")
        header = infile.header.copy()
        outfile = pysam.AlignmentFile(args.output, "wb", header=header)
    except FileNotFoundError:
        print(f"Error: Input BAM file not found at {args.bam}", file=sys.stderr)
        sys.exit(1)
    except ValueError as e:
        print(f"Error opening BAM file {args.bam}. Is it a valid BAM file? Is the index missing?", file=sys.stderr)
        print(f"Pysam error: {e}", file=sys.stderr)
        sys.exit(1)


    # --- 3. Initialize counters and structures for novelty ---
    novel_gene_counter = 1
    novel_transcript_counter = 1
    # This tree tracks the genomic loci of newly discovered multi-exon novel genes
    novel_gene_loci = defaultdict(IntervalTree)

    # --- 4. Process each read in the BAM file ---
    read_count = 0
    annotated_count = 0
    for read in infile:
        read_count += 1
        if read_count % 100000 == 0:
            print(f"[{time.ctime()}] Processed {read_count} reads...")

        # Skip unmapped reads or secondary alignments
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            outfile.write(read)
            continue

        read_chrom = read.reference_name
        read_start = read.reference_start
        read_end = read.reference_end
        read_strand = get_read_strand(read)
        
        # get_blocks() gives the exon structure for spliced reads
        read_blocks = read.get_blocks()
        
        # --- Find candidate gene overlaps using the interval tree ---
        candidate_genes = gene_interval_trees[read_chrom].overlap(read_start, read_end)
        
        # Filter candidates by strand
        stranded_candidates = {
            interval.data for interval in candidate_genes 
            if gene_annotations[interval.data]['strand'] == read_strand
        }

        best_gene_id = None
        assigned_tx_id = None

        if not stranded_candidates:
            # --- NOVEL GENE LOGIC ---
            # This read doesn't overlap any known gene on the same strand.
            # Check if it overlaps a previously defined novel gene locus.
            novel_locus_overlap = novel_gene_loci[read_chrom].overlap(read_start, read_end)
            
            overlapping_novel_genes = {
                interval.data for interval in novel_locus_overlap
                if interval.data[1] == read_strand # data is (novel_id, strand)
            }

            if overlapping_novel_genes:
                # Assign to the first overlapping novel gene found
                best_gene_id = overlapping_novel_genes.pop()[0]
            else:
                # This is a new novel gene locus
                best_gene_id = f"NOVELG{novel_gene_counter:04d}"
                novel_gene_counter += 1
                # Store the novel gene with its strand
                novel_gene_loci[read_chrom].add(Interval(read_start, read_end, (best_gene_id, read_strand)))
            
            # Every novel gene read also gets a novel transcript ID
            assigned_tx_id = f"NOVELT{novel_transcript_counter:012d}"
            novel_transcript_counter += 1

        else:
            # --- KNOWN GENE AND TRANSCRIPT LOGIC ---
            # This read overlaps at least one known gene.
            
            # AMBIGUITY RESOLUTION: Find the gene with the greatest exonic overlap.
            max_overlap = -1
            
            for gene_id in stranded_candidates:
                overlap = calculate_exonic_overlap(read_blocks, gene_annotations[gene_id]['exons'])
                if overlap > max_overlap:
                    max_overlap = overlap
                    best_gene_id = gene_id
            
            if best_gene_id is None:
                # This can happen if the read overlaps the intron of a single-exon gene
                # but has no exonic overlap. Treat as novel.
                # This logic path is rare but necessary for correctness.
                best_gene_id = f"NOVELG{novel_gene_counter:04d}"
                novel_gene_counter += 1
                assigned_tx_id = f"NOVELT{novel_transcript_counter:012d}"
                novel_transcript_counter += 1
            else:
                # We have a best gene. Now check for a perfect transcript match.
                read_exon_tuple = tuple(sorted(read_blocks))
                
                found_match = False
                for tx_id in gene_to_transcripts.get(best_gene_id, []):
                    if transcript_annotations[tx_id]['exons'] == read_exon_tuple:
                        assigned_tx_id = tx_id
                        found_match = True
                        break
                
                if not found_match:
                    # NOVEL TRANSCRIPT, KNOWN GENE
                    assigned_tx_id = f"NOVELT{novel_transcript_counter:012d}"
                    novel_transcript_counter += 1

        # --- Set tags and write the read to the output file ---
        if best_gene_id and assigned_tx_id:
            read.set_tag('gx', best_gene_id, 'Z')
            read.set_tag('tx', assigned_tx_id, 'Z')
            annotated_count += 1
        
        outfile.write(read)

    # --- 5. Clean up and finalize ---
    print(f"[{time.ctime()}] Finished processing {read_count} reads.")
    print(f"[{time.ctime()}] Annotated {annotated_count} reads with gx/tx tags.")
    infile.close()
    outfile.close()

    # --- 6. Sort and index the output BAM file for compatibility ---
    print(f"[{time.ctime()}] Sorting and indexing output BAM file...")
    temp_sorted_path = args.output + ".sorted"
    pysam.sort("-o", temp_sorted_path, args.output)
    os.rename(temp_sorted_path, args.output)
    pysam.index(args.output)
    print(f"[{time.ctime()}] Script finished. Output at: {args.output}")


if __name__ == "__main__":
    main()
