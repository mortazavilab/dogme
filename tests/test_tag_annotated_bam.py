import gzip
import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "tag_annotated_bam.py"
MODULE_SPEC = importlib.util.spec_from_file_location("tag_annotated_bam", SCRIPT)
tag_annotated_bam = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(tag_annotated_bam)

pysam = pytest.importorskip("pysam")


def write_bam(path, read_ids):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    with pysam.AlignmentFile(path, "wb", header=header) as bam:
        for position, read_id in enumerate(read_ids):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = read_id
            read.query_sequence = "A" * 10
            read.flag = 0
            read.reference_id = 0
            read.reference_start = position * 10
            read.mapping_quality = 60
            read.cigar = ((0, 10),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 10)
            read.set_tag("GX", "gene1")
            read.set_tag("TX", "transcript1")
            bam.write(read)


def write_fastq(path, records):
    with gzip.open(path, "wt") as output:
        for read_id, sequence in records:
            output.write(f"@{read_id}\n{sequence}\n+\n{'I' * len(sequence)}\n")


def test_adds_splitcode_tags_and_reports_reads_not_emitted_by_splitcode(tmp_path):
    bam = tmp_path / "input.bam"
    umi = tmp_path / "umi.fastq.gz"
    barcode = tmp_path / "barcode.fastq.gz"
    output = tmp_path / "tagged.bam"
    qc = tmp_path / "tagged.qc.tsv"
    write_bam(bam, ["read1", "read2", "read3"])
    write_fastq(umi, [("read1", "AACCGGTTAA"), ("read2", "TTGGCCAATT")])
    write_fastq(barcode, [("read1", "ACGTACGTACGTACGTACGTACGT"), ("read2", "TGCATGCATGCATGCATGCATGCA")])

    metrics = tag_annotated_bam.tag_bam(bam, umi, barcode, output, qc)

    with pysam.AlignmentFile(output, "rb") as tagged:
        reads = {read.query_name: read for read in tagged.fetch(until_eof=True)}
    assert reads["read1"].get_tag("CB") == "ACGTACGTACGTACGTACGTACGT"
    assert reads["read1"].get_tag("UB") == "AACCGGTTAA"
    assert reads["read1"].get_tag("GX") == "gene1"
    assert reads["read1"].get_tag("TX") == "transcript1"
    assert reads["read2"].get_tag("CB") == "TGCATGCATGCATGCATGCATGCA"
    assert not reads["read3"].has_tag("CB")
    assert metrics == {
        "splitcode_reads": 2,
        "bam_alignments": 3,
        "tagged_bam_alignments": 2,
        "untagged_bam_alignments": 1,
        "splitcode_reads_absent_from_annotated_bam": 0,
    }
    assert "untagged_bam_alignments\t1" in qc.read_text()


def test_rejects_misaligned_splitcode_fastq_records(tmp_path):
    bam = tmp_path / "input.bam"
    umi = tmp_path / "umi.fastq.gz"
    barcode = tmp_path / "barcode.fastq.gz"
    write_bam(bam, ["read1"])
    write_fastq(umi, [("read1", "AACCGGTTAA")])
    write_fastq(barcode, [("read2", "ACGTACGTACGTACGTACGTACGT")])

    with pytest.raises(ValueError, match="out of sync"):
        tag_annotated_bam.tag_bam(
            bam, umi, barcode, tmp_path / "tagged.bam", tmp_path / "qc.tsv"
        )


def test_reports_splitcode_read_ids_missing_from_annotated_bam(tmp_path):
    bam = tmp_path / "input.bam"
    umi = tmp_path / "umi.fastq.gz"
    barcode = tmp_path / "barcode.fastq.gz"
    write_bam(bam, ["read1"])
    write_fastq(umi, [("not_in_bam", "AACCGGTTAA")])
    write_fastq(barcode, [("not_in_bam", "ACGTACGTACGTACGTACGTACGT")])

    qc = tmp_path / "qc.tsv"
    metrics = tag_annotated_bam.tag_bam(
        bam, umi, barcode, tmp_path / "tagged.bam", qc
    )

    assert metrics["tagged_bam_alignments"] == 0
    assert metrics["splitcode_reads_absent_from_annotated_bam"] == 1
    assert "splitcode_reads_absent_from_annotated_bam\t1" in qc.read_text()


def test_rejects_duplicate_splitcode_read_ids(tmp_path):
    bam = tmp_path / "input.bam"
    umi = tmp_path / "umi.fastq.gz"
    barcode = tmp_path / "barcode.fastq.gz"
    write_bam(bam, ["read1"])
    write_fastq(umi, [("read1", "AACCGGTTAA"), ("read1", "TTGGCCAATT")])
    write_fastq(barcode, [("read1", "ACGTACGTACGTACGTACGTACGT"), ("read1", "TGCATGCATGCATGCATGCATGCA")])

    with pytest.raises(ValueError, match="Duplicate read ID"):
        tag_annotated_bam.tag_bam(
            bam, umi, barcode, tmp_path / "tagged.bam", tmp_path / "qc.tsv"
        )
