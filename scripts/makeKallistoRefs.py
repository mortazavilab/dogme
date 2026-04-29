#!/usr/bin/env python3
import sys
import argparse
import re
from collections import defaultdict

def parse_args():
    parser = argparse.ArgumentParser(description="Generate Kallisto Total index files from FASTA and GTF.")
    parser.add_argument('--name', required=True, help="Prefix for output files (e.g., species1)")
    parser.add_argument('--fasta', required=True, help="Reference genome FASTA file")
    parser.add_argument('--gtf', required=True, help="Annotation GTF file")
    return parser.parse_args()

def get_reverse_complement(seq):
    # Handles standard and ambiguous nucleotides
    trans = str.maketrans('ATGCUatgcuNn', 'TACGAtacgaNn')
    return seq.translate(trans)[::-1]

def parse_gtf(gtf_path):
    """
    Parses the GTF to extract exon coordinates grouped by chromosome and transcript_id.
    Also builds the transcript-to-gene mapping.
    """
    chrom_data = defaultdict(lambda: defaultdict(list))
    t2g = {}
    
    # Regex to safely extract IDs from the attributes column
    transcript_re = re.compile(r'transcript_id "([^"]+)"')
    gene_re = re.compile(r'gene_id "([^"]+)"')
    
    print(f"Parsing GTF: {gtf_path}...", file=sys.stderr)
    with open(gtf_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            # We only care about exon features for transcript building
            if len(parts) < 9 or parts[2] != 'exon':
                continue
            
            chrom = parts[0]
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            attributes = parts[8]
            
            t_match = transcript_re.search(attributes)
            g_match = gene_re.search(attributes)
            
            if not t_match or not g_match:
                continue
                
            t_id = t_match.group(1)
            g_id = g_match.group(1)
            
            t2g[t_id] = g_id
            chrom_data[chrom][t_id].append({
                'start': start,
                'end': end,
                'strand': strand
            })
            
    return chrom_data, t2g

def process_chromosome(chrom, seq_list, chrom_transcripts, cdna_out, intron_out):
    """
    Stitches exons into cDNA and calculates intronic spaces for a single chromosome.
    """
    if not chrom_transcripts:
        return
        
    chrom_seq = "".join(seq_list)
    seq_length = len(chrom_seq)
    
    for t_id, exons in chrom_transcripts.items():
        # Sort exons strictly by genomic start coordinate (ascending)
        exons.sort(key=lambda x: x['start'])
        strand = exons[0]['strand']
        
        cdna_seqs = []
        # Extract and build cDNA
        for ex in exons:
            # GTF is 1-based inclusive. Python strings are 0-based exclusive.
            # Example: GTF 10 to 15 -> Python slice [9:15]
            ex_seq = chrom_seq[ex['start']-1 : ex['end']]
            cdna_seqs.append(ex_seq)
        
        full_cdna = "".join(cdna_seqs)
        
        # If on the negative strand, the biological RNA is the reverse complement
        if strand == '-':
            full_cdna = get_reverse_complement(full_cdna)
            
        # Write cDNA wrapped to 80 characters
        cdna_out.write(f">{t_id}\n")
        for i in range(0, len(full_cdna), 80):
            cdna_out.write(full_cdna[i:i+80] + "\n")
            
        # Extract and build Introns (the spaces between sorted exons)
        for i in range(len(exons) - 1):
            intron_start = exons[i]['end']  # End of current exon
            intron_end = exons[i+1]['start'] - 1  # Start of next exon - 1
            
            if intron_end > intron_start:
                intron_seq = chrom_seq[intron_start : intron_end]
                if strand == '-':
                    intron_seq = get_reverse_complement(intron_seq)
                    
                # Ensure unique FASTA headers for Kallisto
                intron_out.write(f">{t_id}_intron_{i+1}\n")
                for j in range(0, len(intron_seq), 80):
                    intron_out.write(intron_seq[j:j+80] + "\n")

def main():
    args = parse_args()
    chrom_data, t2g = parse_gtf(args.gtf)
    
    # Write t2g mapping file
    t2g_file = f"{args.name}.t2g"
    print(f"Writing {t2g_file}...", file=sys.stderr)
    with open(t2g_file, 'w') as f:
        for t_id, g_id in t2g.items():
            f.write(f"{t_id}\t{g_id}\n")
            
    # Process FASTA file chromosome by chromosome
    cdna_file = f"{args.name}.cdna.fa"
    intron_file = f"{args.name}.introns.fa"
    
    print(f"Streaming FASTA: {args.fasta}...", file=sys.stderr)
    with open(args.fasta, 'r') as fasta_in, \
         open(cdna_file, 'w') as cdna_out, \
         open(intron_file, 'w') as intron_out:
         
        current_chrom = None
        current_seq = []
        
        for line in fasta_in:
            line = line.strip()
            if line.startswith('>'):
                # We hit a new chromosome header; process the previous one if it exists
                if current_chrom is not None:
                    if current_chrom in chrom_data:
                        print(f"  Extracting sequences for {current_chrom}...", file=sys.stderr)
                        process_chromosome(current_chrom, current_seq, chrom_data[current_chrom], cdna_out, intron_out)
                    else:
                        print(f"  Skipping {current_chrom} (not found in GTF)...", file=sys.stderr)
                
                # Start collecting the new chromosome (extracting just the ID up to the first space)
                current_chrom = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
                
        # Don't forget to process the final chromosome in the file
        if current_chrom is not None and current_chrom in chrom_data:
            print(f"  Extracting sequences for {current_chrom}...", file=sys.stderr)
            process_chromosome(current_chrom, current_seq, chrom_data[current_chrom], cdna_out, intron_out)
            
    print("Done! Files are ready for Kallisto indexing.", file=sys.stderr)

if __name__ == '__main__':
    main()
