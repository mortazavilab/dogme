import sys

def parse_gtf(gtf_file):
    with open(gtf_file) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if fields[2] != "exon":
                continue
            chrom = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            transcript_id = None
            for attr in fields[8].split(';'):
                if "transcript_id" in attr:
                    transcript_id = attr.split('"')[1]
            yield chrom, start, end, strand, transcript_id

def get_junctions(gtf_file):
    transcripts = {}
    for chrom, start, end, strand, tid in parse_gtf(gtf_file):
        if tid not in transcripts:
            transcripts[tid] = []
        transcripts[tid].append((chrom, start, end, strand))
    for tid, exons in transcripts.items():
        exons.sort(key=lambda x: x[1])
        for i in range(len(exons) - 1):
            chrom, _, exon_end, strand = exons[i]
            _, next_exon_start, _, _ = exons[i+1]
            # BED is 0-based, half-open
            bed_start = exon_end
            bed_end = next_exon_start - 1
            if bed_end > bed_start:
                print(f"{chrom}\t{bed_start}\t{bed_end}\t{tid}\t0\t{strand}")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python gtf_to_junction_bed.py <input.gtf>", file=sys.stderr)
        sys.exit(1)
    get_junctions(sys.argv[1])