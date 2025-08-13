#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Long-Read RNA-Seq BAM Annotation Script

This script annotates long-read RNA-sequencing alignments in a BAM file
based on a provided GTF/GFF3 gene annotation. It classifies each read into
novelty categories (Known, ISM, NIC, NNC, Antisense, Intergenic) based on
its splice junction compatibility with the known annotation.

The script is optimized for performance using multiprocessing and produces
an annotated BAM file with custom tags. It includes a tolerance for "fuzzy"
junctions common in long-read alignments.
"""

import argparse
import sys
import time
import re
import os
import subprocess
from collections import defaultdict, Counter
from multiprocessing import Pool, Manager
from multiprocessing.managers import SyncManager
import pysam

__version__ = "1.1.0" 
__author__ = "Elnaz A., Gemini, Claude, Ali M."

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Logging ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def log_message(message, level="INFO"):
    """Prints a formatted log message to stderr."""
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    sys.stderr.write(f"[{timestamp}][{level}] {message}\n")
    sys.stderr.flush()

# ~~~~~~~~~~~~~~~~~~~~~~ GTF/GFF3 Parsing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def parse_gtf(gtf_file):
    """
    Parses a GTF/GFF3 file and builds annotation data structures.
    """
    log_message(f"Parsing annotation file: {gtf_file}")
    genes = defaultdict(dict)
    all_donors = defaultdict(lambda: defaultdict(set))
    all_acceptors = defaultdict(lambda: defaultdict(set))
    gene_id_to_name = {}
    gene_id_re = re.compile(r'gene_id[ =]"?([^";]+)"?')
    transcript_id_re = re.compile(r'transcript_id[ =]"?([^";]+)"?')
    gene_name_re = re.compile(r'gene_name[ =]"?([^";]+)"?')

    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            if feature_type not in ('exon', 'transcript', 'gene'):
                continue

            chrom, start, end, strand = fields[0], int(fields[3]), int(fields[4]), fields[6]
            attributes = fields[8]

            gene_id_match = gene_id_re.search(attributes)
            transcript_id_match = transcript_id_re.search(attributes)
            gene_name_match = gene_name_re.search(attributes)

            if not gene_id_match:
                continue
            gene_id = gene_id_match.group(1).split('.')[0]
            if gene_name_match:
                gene_id_to_name[gene_id] = gene_name_match.group(1)

            if gene_id not in genes[chrom]:
                genes[chrom][gene_id] = {
                    'strand': strand,
                    'start': float('inf'),
                    'end': float('-inf'),
                    'transcripts': defaultdict(lambda: {'exons': []}),
                    'donors': set(),
                    'acceptors': set()
                }
            genes[chrom][gene_id]['start'] = min(genes[chrom][gene_id]['start'], start)
            genes[chrom][gene_id]['end'] = max(genes[chrom][gene_id]['end'], end)

            if feature_type == 'exon' and transcript_id_match:
                transcript_id = transcript_id_match.group(1).split('.')[0]
                genes[chrom][gene_id]['transcripts'][transcript_id]['exons'].append((start, end))

    for chrom in genes:
        for gene_id in genes[chrom]:
            gene = genes[chrom][gene_id]
            strand = gene['strand']
            for transcript_id in gene['transcripts']:
                gene['transcripts'][transcript_id]['exons'].sort()
                exons = gene['transcripts'][transcript_id]['exons']
                for i in range(len(exons) - 1):
                    # For '-' strand, transcription flows from higher to lower coordinates.
                    # Exon[i+1] is transcribed before exon[i].
                    # Donor is at start of exon[i+1], Acceptor is at end of exon[i].
                    if strand == '+':
                        donor, acceptor = exons[i][1], exons[i+1][0]
                        all_donors[chrom][strand].add(donor)
                        all_acceptors[chrom][strand].add(acceptor)
                        gene['donors'].add(donor)
                        gene['acceptors'].add(acceptor)
                    else: # strand == '-'
                        donor, acceptor = exons[i+1][0], exons[i][1]
                        all_donors[chrom][strand].add(donor)
                        all_acceptors[chrom][strand].add(acceptor)
                        gene['donors'].add(donor)
                        gene['acceptors'].add(acceptor)

    log_message("Finished parsing annotation.")
    return genes, all_donors, all_acceptors, gene_id_to_name

# ~~~~~~~~~~~~~~~~~~~~~~ Core Annotation Logic ~~~~~~~~~~~~~~~~~~~~~~~~~~
def get_read_splice_junctions(read, min_intron_len=10):
    """
    Extracts splice junctions from a BAM alignment record.
    Coordinates are 1-based, inclusive, from 0-based pysam blocks.
    """
    junctions = []
    blocks = read.get_blocks()
    if len(blocks) < 2:
        return []

    for i in range(len(blocks) - 1):
        donor = blocks[i][1]
        acceptor = blocks[i+1][0] + 1
        intron_length = blocks[i+1][0] - blocks[i][1]
        if intron_length >= min_intron_len:
            junctions.append((donor, acceptor))
    return junctions

def classify_read(read, annotation, min_intron_len, junc_tolerance, debug_gene=None, debug_list=None): # <<< CHANGED: Added junc_tolerance
    """
    Classifies a single read based on its splice junctions against the annotation,
    allowing for a configurable tolerance at junction sites.
    """
    if read.is_unmapped or read.is_secondary or read.is_supplementary:
        return None

    genes, all_donors, all_acceptors = annotation
    read_chrom = read.reference_name
    read_start = read.reference_start + 1
    read_end = read.reference_end
    read_strand = '-' if read.is_reverse else '+'
    read_junctions = get_read_splice_junctions(read, min_intron_len)
    overlapping_genes_same_strand = []
    overlapping_genes_opp_strand = []

    if read_chrom in genes:
        for gene_id, gene_data in genes[read_chrom].items():
            if read_start < gene_data['end'] and read_end > gene_data['start']:
                if gene_data['strand'] == read_strand:
                    overlapping_genes_same_strand.append((gene_id, gene_data))
                else:
                    overlapping_genes_opp_strand.append((gene_id, gene_data))

    result = None
    if overlapping_genes_same_strand:
        potential_matches = []
        for gene_id, gene_data in overlapping_genes_same_strand:
            for transcript_id, trans_data in gene_data['transcripts'].items():
                known_exons = trans_data['exons']
                known_junctions = [(known_exons[i][1], known_exons[i+1][0]) for i in range(len(known_exons) - 1)]

                # <<< CHANGED: Logic replaced with tolerant matching
                # Check for a full splice match within tolerance
                is_full_match = False
                if len(read_junctions) == len(known_junctions):
                    is_full_match = all(
                        abs(r_j[0] - k_j[0]) <= junc_tolerance and abs(r_j[1] - k_j[1]) <= junc_tolerance
                        for r_j, k_j in zip(read_junctions, known_junctions)
                    )
                
                if is_full_match:
                    potential_matches.append({"is_full_match": True, "gene_id": gene_id, "transcript_id": transcript_id})
                    continue # Found best possible match for this transcript

                # Check for an incomplete splice match (subset) within tolerance
                is_subset_match = False
                if read_junctions:
                    matched_known_indices = set()
                    num_read_juncs_matched = 0
                    for r_j in read_junctions:
                        for i, k_j in enumerate(known_junctions):
                            if i not in matched_known_indices and (abs(r_j[0] - k_j[0]) <= junc_tolerance and abs(r_j[1] - k_j[1]) <= junc_tolerance):
                                matched_known_indices.add(i)
                                num_read_juncs_matched += 1
                                break
                    if num_read_juncs_matched == len(read_junctions):
                        is_subset_match = True

                if is_subset_match:
                    potential_matches.append({"is_full_match": False, "gene_id": gene_id, "transcript_id": transcript_id})

        full_matches = [m for m in potential_matches if m['is_full_match']]
        ism_matches = [m for m in potential_matches if not m['is_full_match']]

        if full_matches:
            best_match = full_matches[0]
            classification = "Known"
            result = (classification, best_match['gene_id'], best_match['transcript_id'])
        elif ism_matches:
            best_match = ism_matches[0]
            classification = "ISM"
            result = (classification, best_match['gene_id'], best_match['transcript_id'])
        else: # Fall through to NIC/NNC logic
            best_gene_id = overlapping_genes_same_strand[0][0] # Default to first overlapping
            if not read_junctions:
                classification = "NIC"
            else:
                is_nic = True
                strand_donors = all_donors.get(read_chrom, {}).get(read_strand, set())
                strand_acceptors = all_acceptors.get(read_chrom, {}).get(read_strand, set())
                for junc in read_junctions:
                    # <<< CHANGED: Tolerant check for NIC/NNC
                    read_d, read_a = (junc[0], junc[1])
                    donor = read_d if read_strand == '+' else read_a
                    acceptor = read_a if read_strand == '+' else read_d

                    donor_found = any((donor + i) in strand_donors for i in range(-junc_tolerance, junc_tolerance + 1))
                    acceptor_found = any((acceptor + i) in strand_acceptors for i in range(-junc_tolerance, junc_tolerance + 1))
                    
                    if not (donor_found and acceptor_found):
                        is_nic = False
                        break
                classification = "NIC" if is_nic else "NNC"
            result = (classification, best_gene_id, "NOVEL")

    elif overlapping_genes_opp_strand:
        gene_id = overlapping_genes_opp_strand[0][0]
        result = ("Antisense", gene_id, "ANTISENSE")
    else:
        result = ("Intergenic", "NOVELG", "NOVELT")

    classification, final_gene_id, final_transcript_id = result

    if debug_gene and final_gene_id == debug_gene and debug_list is not None:
        log_prefix = f"Read:{read.query_name}\tJunctions:{read_junctions}\t->\tClass:{classification}\tGene:{final_gene_id}\tTranscript:{final_transcript_id}"
        debug_list.append(log_prefix)

    read.set_tag("TT", classification)
    read.set_tag("GX", final_gene_id)
    read.set_tag("TX", final_transcript_id)
    return classification, final_gene_id, final_transcript_id, read

# ~~~~~~~~~~~~~~~~~~~~~~ Parallel Processing & I/O ~~~~~~~~~~~~~~~~~~~~~~
# <<< CHANGED: Added junc_tolerance to worker initialization
def worker_init(annotation_data, header_dict, min_len, tolerance, debug_gene_id=None, debug_novel_type=None, debug_gene_list=None, debug_novel_list=None, debug_gene_info=None):
    global worker_annotation, worker_header, worker_min_intron_len, worker_junc_tolerance, worker_debug_gene, worker_debug_novel_type, worker_debug_gene_list, worker_debug_novel_list, worker_debug_gene_info
    worker_annotation = annotation_data
    worker_header = pysam.AlignmentHeader.from_dict(header_dict)
    worker_min_intron_len = min_len
    worker_junc_tolerance = tolerance # <<< CHANGED
    worker_debug_gene = debug_gene_id
    worker_debug_novel_type = debug_novel_type
    worker_debug_gene_list = debug_gene_list
    worker_debug_novel_list = debug_novel_list
    worker_debug_gene_info = debug_gene_info

def process_chunk(read_strings):
    results = []
    for read_string in read_strings:
        read = pysam.AlignedSegment.fromstring(read_string, worker_header)
        # <<< CHANGED: Pass junc_tolerance to classify_read
        result = classify_read(read, worker_annotation, worker_min_intron_len, worker_junc_tolerance, worker_debug_gene, worker_debug_gene_list)
        if result:
            classification, gene_id, transcript_id, annotated_read = result
            results.append((classification, gene_id, transcript_id, annotated_read.to_string()))
            if worker_debug_novel_type and classification == worker_debug_novel_type:
                junctions = get_read_splice_junctions(read, worker_min_intron_len)
                debug_line = f"Read:{read.query_name}\tGene:{gene_id}\tTranscript:{transcript_id}\tJunctions:{junctions}"
                worker_debug_novel_list.append(debug_line)
    return results

def read_bam_in_chunks(bam_file, chunk_size, max_reads=None, region=None):
    read_iterator = bam_file.fetch(**region) if region else bam_file.fetch(until_eof=True)
    read_chunk = []
    reads_yielded = 0
    for read in read_iterator:
        if max_reads is not None and reads_yielded >= max_reads:
            break
        read_chunk.append(read.to_string())
        reads_yielded += 1
        if len(read_chunk) >= chunk_size:
            yield read_chunk
            read_chunk = []
    if read_chunk:
        yield read_chunk

# ~~~~~~~~~~~~~~~~~~~~~~ Output Generation (No changes needed here) ~~~~~
def generate_qc_report(abundance_counter, gene_id_to_name, qc_filename):
    # This function is correct as is.
    log_message(f"Generating QC report: {qc_filename}")
    gene_qc_data = defaultdict(lambda: defaultdict(lambda: {'models': set(), 'reads': 0}))
    for (gene_id, transcript_id, classification), count in abundance_counter.items():
        if count > 0:
            gene_qc_data[gene_id][classification]['models'].add(transcript_id)
            gene_qc_data[gene_id][classification]['reads'] += count
    with open(qc_filename, 'w') as f:
        header = [
            "gene_id", "gene_name", "known_models", "known_reads", "ism_models", "ism_reads",
            "nic_models", "nic_reads", "nnc_models", "nnc_reads", "antisense_models", "antisense_reads",
            "intergenic_models", "intergenic_reads"
        ]
        f.write(",".join(header) + "\n")
        for gene_id in sorted(gene_qc_data.keys()):
            gene_name = gene_id_to_name.get(gene_id, 'NovelGene')
            def get_counts(classification):
                data = gene_qc_data[gene_id].get(classification, {})
                return len(data.get('models', set())), data.get('reads', 0)

            known_models, known_reads = get_counts('Known')
            ism_models, ism_reads = get_counts('ISM')
            nic_models, nic_reads = get_counts('NIC')
            nnc_models, nnc_reads = get_counts('NNC')
            antisense_models, antisense_reads = get_counts('Antisense')
            intergenic_models, intergenic_reads = get_counts('Intergenic')
            f.write(f"{gene_id},{gene_name},{known_models},{known_reads},{ism_models},{ism_reads},"
                    f"{nic_models},{nic_reads},{nnc_models},{nnc_reads},{antisense_models},{antisense_reads},"
                    f"{intergenic_models},{intergenic_reads}\n")

def sort_and_index_bam(unsorted_bam_path, sorted_bam_path):
    # This function is correct as is.
    log_message(f"Sorting BAM file: {unsorted_bam_path}")
    try:
        subprocess.run(['samtools', 'sort', '-@', '4', '-o', sorted_bam_path, unsorted_bam_path], check=True, stderr=subprocess.DEVNULL)
        log_message(f"Indexing BAM file: {sorted_bam_path}")
        subprocess.run(['samtools', 'index', sorted_bam_path], check=True)
        os.remove(unsorted_bam_path)
        log_message(f"Removed temporary unsorted file: {unsorted_bam_path}")
    except FileNotFoundError:
        log_message("ERROR: samtools not found. Please ensure it is installed and in your PATH.", level="ERROR")
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        log_message(f"ERROR: samtools command failed: {e}", level="ERROR")
        sys.exit(1)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~ Main Orchestrator ~~~~~~~~~~~~~~~~~~~~~~~~~~~
def main():
    parser = argparse.ArgumentParser(description="Annotate long-read RNA-seq BAM files.", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("--bam", "-b", required=True, help="Input BAM file.")
    parser.add_argument("--gtf", "-g", required=True, help="Gene annotation file.")
    parser.add_argument("--out", "-o", required=True, help="Output file prefix.")
    parser.add_argument("--threads", "-t", type=int, default=1, help="Number of threads.")
    parser.add_argument("--talon", action="store_true", help="Generate TALON-compatible output.")
    parser.add_argument("--chunk_size", type=int, default=5000, help="Reads per chunk.")
    parser.add_argument("--num_reads", type=int, default=None, help="Process first N reads.")
    parser.add_argument("--min_intron_len", type=int, default=10, help="Min intron length.")
    # <<< CHANGED: Added junction tolerance argument
    parser.add_argument("--junc_tolerance", type=int, default=2, help="Tolerance for matching splice junctions (e.g., 2 for +/- 2 bp).")
    parser.add_argument("--debug_gene", type=str, default=None, help="Debug a single gene ID, writing output to a file.")
    parser.add_argument("--debug_novel_type", type=str, default=None, help="Output details for a specific novel class (e.g., NNC, NIC).")
    args = parser.parse_args()

    start_time = time.time()
    if args.debug_gene:
        args.threads = 1
        log_message(f"DEBUG MODE: Focusing on gene '{args.debug_gene}'. Forcing threads to 1.", level="DEBUG")

    log_message(f"Starting annotation process with {args.threads} threads.")
    log_message(f"Junction tolerance set to: +/- {args.junc_tolerance} bp.") # <<< CHANGED
    if args.num_reads:
        log_message(f"Processing only the first {args.num_reads:,} reads.")
    
    genes, all_donors, all_acceptors, gene_id_to_name = parse_gtf(args.gtf)
    debug_region = None
    debug_gene_info = None

    if args.debug_gene:
        # Debug logic is fine
        found = False
        for chrom, gene_dict in genes.items():
            if args.debug_gene in gene_dict:
                gene_info = gene_dict[args.debug_gene]
                debug_gene_info = gene_info
                debug_region = {'contig': chrom, 'start': gene_info['start'] - 1, 'end': gene_info['end']}
                log_message(f"Found debug gene {args.debug_gene} at {chrom}:{gene_info['start']}-{gene_info['end']}", level="DEBUG")
                found = True
                break
        if not found:
            log_message(f"ERROR: Debug gene '{args.debug_gene}' not found in GTF.", level="ERROR")
            sys.exit(1)

    input_bam = pysam.AlignmentFile(args.bam, "rb")
    if args.out.lower().endswith('.bam'):
        # If user provides a full .bam filename, use it
        sorted_bam_file = args.out
        # And derive the prefix from it for other files
        output_prefix = os.path.splitext(sorted_bam_file)[0] 
    else:
        # Otherwise, assume a prefix was given and build the name
        output_prefix = args.out
        sorted_bam_file = f"{output_prefix}.annotated.bam"

    # The other files are now based on the derived prefix
    unsorted_bam_file = f"{output_prefix}.annot-unsorted.bam"
    qc_report_file = f"{output_prefix}_qc_summary.csv"

    log_message(f"Temporary unsorted BAM will be written to: {unsorted_bam_file}")
    output_bam = pysam.AlignmentFile(unsorted_bam_file, "wb", template=input_bam)
    header_dict = input_bam.header.to_dict()

    manager = Manager()
    debug_gene_log_list = manager.list() if args.debug_gene else None
    debug_novel_log_list = manager.list() if args.debug_novel_type else None
    # <<< CHANGED: Pass junction tolerance to workers
    pool = Pool(processes=args.threads, initializer=worker_init, initargs=(
        (genes, all_donors, all_acceptors), header_dict, args.min_intron_len, args.junc_tolerance,
        args.debug_gene, args.debug_novel_type, debug_gene_log_list, debug_novel_log_list, debug_gene_info
    ))

    abundance_counter = Counter()
    classification_counter = Counter()
    known_genes_found = {}
    known_transcripts_found = {}
    novel_transcripts = {}
    novel_genes = {}
    novel_gene_counter = 0
    novel_transcript_counter = 0

    log_message("Submitting and processing reads...")
    main_process_header = pysam.AlignmentHeader.from_dict(header_dict)
    chunk_generator = read_bam_in_chunks(input_bam, args.chunk_size, args.num_reads, debug_region)
    total_reads_processed = 0
    next_progress_report = 100000

    # Main processing loop has been simplified to not use multiprocessing Counters
    # for better performance and simplicity, aggregation happens at the end.
    all_results = []
    for chunk_results in pool.imap_unordered(process_chunk, chunk_generator):
        all_results.extend(chunk_results)
        total_reads_processed += len(chunk_results)
        if not args.debug_gene and total_reads_processed >= next_progress_report:
            log_message(f"Processed {total_reads_processed:,} reads...")
            next_progress_report += 100000

    pool.close()
    pool.join()
    
    log_message("Aggregating results and writing final BAM...")
    for classification, gene_id, transcript_id, read_string in all_results:
        read = pysam.AlignedSegment.fromstring(read_string, main_process_header)
        read_strand = '-' if read.is_reverse else '+'
        final_gene_id, final_transcript_id = gene_id, transcript_id

        if gene_id == "NOVELG":
            gene_locus = (read.reference_name, read.reference_start, read.reference_end, read_strand)
            if gene_locus not in novel_genes:
                novel_gene_counter += 1
                novel_genes[gene_locus] = f"NOVELG{novel_gene_counter:06d}"
            final_gene_id = novel_genes[gene_locus]
            read.set_tag("GX", final_gene_id)

        if transcript_id in ["NOVEL", "NOVELT"]:
            transcript_locus = (final_gene_id, tuple(get_read_splice_junctions(read, args.min_intron_len)), read_strand)
            if transcript_locus not in novel_transcripts:
                novel_transcript_counter += 1
                new_id = f"NOVELT{novel_transcript_counter:010d}"
                novel_transcripts[transcript_locus] = new_id
            final_transcript_id = novel_transcripts[transcript_locus]
            read.set_tag("TX", final_transcript_id)

        output_bam.write(read)

        if gene_id != "NOVELG":
            known_genes_found[final_gene_id] = True
        if classification in ["Known", "ISM"]:
            known_transcripts_found[final_transcript_id] = True
        
        abundance_key = (final_gene_id, final_transcript_id, classification)
        abundance_counter.update([abundance_key])
        classification_counter.update([classification])

    input_bam.close()
    output_bam.close()
    
    # Debugging and report generation logic is fine, no changes needed
    if args.debug_gene and debug_gene_log_list:
        debug_filename = f"{output_prefix}_debug_gene_{args.debug_gene}.txt"
        log_message(f"Writing gene debug log to {debug_filename}")
        with open(debug_filename, 'w') as f:
            if debug_gene_info:
                 for tx_id, tx_data in sorted(debug_gene_info['transcripts'].items()):
                    exons = tx_data['exons']
                    junctions = [(exons[i][1], exons[i+1][0]) for i in range(len(exons) - 1)]
                    f.write(f"GTF_MODEL\tTranscript:{tx_id}\tJunctions:{junctions}\n")
            f.write("\n--- READ CLASSIFICATIONS ---\n")
            for line in sorted(list(debug_gene_log_list)):
                f.write(line + '\n')

    if args.debug_novel_type and debug_novel_log_list:
        debug_filename = f"{output_prefix}_debug_{args.debug_novel_type}.txt"
        log_message(f"Writing debug log for '{args.debug_novel_type}' to {debug_filename}")
        with open(debug_filename, 'w') as f:
            for line in sorted(list(debug_novel_log_list)):
                f.write(line + '\n')
    
    generate_qc_report(abundance_counter, gene_id_to_name, qc_report_file)
    sort_and_index_bam(unsorted_bam_file, sorted_bam_file)

    log_message("--- Final Statistics ---")
    log_message(f"Known genes detected: {len(known_genes_found)}")
    log_message(f"Known transcripts detected: {len(known_transcripts_found)}")
    log_message(f"Novel genes discovered: {len(novel_genes)}")
    log_message(f"Novel transcripts discovered: {len(novel_transcripts)}")
    log_message("--- Read Classification Summary ---")
    for classification_type in sorted(classification_counter.keys()):
        count = classification_counter[classification_type]
        log_message(f"{classification_type}: {count:,} reads")
    log_message("------------------------")
    end_time = time.time()
    log_message(f"Annotation complete. Final output file: {sorted_bam_file}")
    log_message(f"Total runtime: {end_time - start_time:.2f} seconds.")

if __name__ == "__main__":
    main()