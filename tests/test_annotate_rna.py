import importlib.util
from collections import Counter
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "annotateRNA.py"
module_spec = importlib.util.spec_from_file_location("annotateRNA", SCRIPT)
annotate_rna = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(annotate_rna)


def test_filter_novel_models_excludes_antisense_transcripts_only_when_disabled():
    novel_models = {
        "known": {"id": "known_tx", "strand": "+"},
        "antisense": {"id": "antisense_tx", "strand": "-"},
    }
    abundance_counter = Counter({
        ("gene1", "known_tx", "Known"): 3,
        ("gene2", "antisense_tx", "Antisense"): 2,
    })

    assert annotate_rna.filter_novel_models_for_export(
        novel_models, abundance_counter, include_antisense=False
    ) == {"known": novel_models["known"]}
    assert annotate_rna.filter_novel_models_for_export(
        novel_models, abundance_counter, include_antisense=True
    ) == novel_models


def make_spliced_read(pysam, cigar):
    header = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6"}, "SQ": [{"SN": "chr1", "LN": 200}]}
    )
    read = pysam.AlignedSegment(header)
    read.query_name = "read"
    read.query_sequence = "A" * sum(length for operation, length in cigar if operation == 0)
    read.flag = 0
    read.reference_id = 0
    read.reference_start = 0
    read.cigar = cigar
    return read


def test_short_single_intron_is_classified_as_genomic():
    pysam = pytest.importorskip("pysam")
    read = make_spliced_read(pysam, [(0, 10), (3, 6), (0, 10)])
    annotation = (
        {
            "chr1": {
                "gene1": {
                    "strand": "+",
                    "start": 1,
                    "end": 26,
                    "transcripts": {"tx1": {"exons": [(1, 10), (17, 26)]}},
                }
            }
        },
        {},
        {},
    )

    assert annotate_rna.get_read_splice_junctions(read, 70) == []
    result = annotate_rna.classify_read(read, annotation, 70, 2)
    assert result[:3] == ("Genomic", "gene1", "tx1")


@pytest.mark.parametrize(
    ("donor", "acceptor", "expected_strand"),
    [("GT", "AG", "+"), ("CT", "AC", "-")],
)
def test_splice_site_motifs_correct_strand_and_report_canonicality(
    tmp_path, donor, acceptor, expected_strand
):
    pysam = pytest.importorskip("pysam")
    sequence = list("A" * 200)
    sequence[20:22] = donor
    sequence[38:40] = acceptor
    fasta_path = tmp_path / "reference.fa"
    fasta_path.write_text(">chr1\n" + "".join(sequence) + "\n")
    pysam.faidx(str(fasta_path))
    fasta = pysam.FastaFile(str(fasta_path))
    read = make_spliced_read(pysam, [(0, 20), (3, 20), (0, 20)])
    read.flag = 16

    assert annotate_rna.infer_splice_strand(read, fasta) == expected_strand
    assert annotate_rna.get_effective_strand(read, fasta) == expected_strand
    assert annotate_rna.get_splice_site_canonicality(
        fasta, "chr1", [(1, 20), (41, 60)]
    ) == "canonical"
    reference_gtf = tmp_path / "reference.gtf"
    flagged_gtf = tmp_path / "flagged.gtf"
    reference_gtf.write_text(
        'chr1\ttest\ttranscript\t1\t60\t.\t+\t.\tgene_id "gene1"; '
        'transcript_id "tx1";\n'
    )
    genes = {
        "chr1": {
            "gene1": {
                "transcripts": {"tx1": {"exons": [(1, 20), (41, 60)]}}
            }
        }
    }
    transcript_info = {"tx1": {}}
    annotate_rna.add_splice_canonicality_to_gtf(
        reference_gtf, flagged_gtf, genes, transcript_info, fasta
    )
    assert 'splice_site_canonical "canonical"' in flagged_gtf.read_text()
    assert transcript_info["tx1"]["splice_site_canonical"] == "canonical"
    fasta.close()


@pytest.mark.parametrize(
    ("repeat_interval", "expected"),
    [("chr1\t18\t42\n", True), ("chr1\t20\t39\n", False)],
)
def test_repeat_filter_requires_complete_intron_containment(
    tmp_path, repeat_interval, expected
):
    pysam = pytest.importorskip("pysam")
    repeat_bed = tmp_path / "simpleRepeat.bed"
    repeat_bed.write_text(repeat_interval)
    repeat_index = annotate_rna.load_simple_repeats(str(repeat_bed))
    read = make_spliced_read(pysam, [(0, 20), (3, 20), (0, 20)])

    assert annotate_rna.junction_fully_contained_in_repeat(
        "chr1", (20, 41), repeat_index
    ) is expected
    assert annotate_rna.get_read_splice_junctions(
        read, 1, repeat_index
    ) == ([] if expected else [(20, 41)])