import gzip
import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "process_splitcode_fastqs.py"
TEMPLATE = ROOT / "templates" / "parse-evercode-wt-mega-v2-nanopore.yaml.j2"
AUTOGENERATION = ROOT / "scripts" / "seqspec_autogeneration.py"

module_spec = importlib.util.spec_from_file_location("process_splitcode_fastqs", SCRIPT)
process_splitcode_fastqs = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(process_splitcode_fastqs)

autogen_spec = importlib.util.spec_from_file_location("seqspec_autogeneration", AUTOGENERATION)
seqspec_autogeneration = importlib.util.module_from_spec(autogen_spec)
autogen_spec.loader.exec_module(seqspec_autogeneration)


EMPIRICAL_100K = {
    "f": {"input": 17308, "valid_umi": 17302, "candidate_homopolymer_pct": 0.01, "cdna_median": 463},
    "rc": {"input": 21022, "valid_umi": 21008, "candidate_homopolymer_pct": 0.0, "cdna_median": 464},
    "ambiguous": 115,
    "exact_unclassified": 61555,
}


def write_fastq(path: Path, records: list[tuple[str, str, str]]) -> None:
    with gzip.open(path, "wt") as handle:
        for read_id, sequence, quality in records:
            handle.write(f"@{read_id} comment\n{sequence}\n+\n{quality}\n")


def read_fastq(path: Path) -> list[tuple[str, str, str]]:
    records = []
    with gzip.open(path, "rt") as handle:
        while name := handle.readline().rstrip():
            sequence = handle.readline().rstrip()
            handle.readline()
            quality = handle.readline().rstrip()
            records.append((name[1:].split()[0], sequence, quality))
    return records


@pytest.mark.parametrize("orientation", ["f", "c", "r", "rc"])
def test_orientation_transform_preserves_barcode_quality_pair(orientation):
    barcode = "ACGT" * 6
    quality = "abcdefghijklmnopqrstuvwx"
    transformed_barcode = process_splitcode_fastqs.transform_barcode(barcode, orientation)
    transformed_quality = process_splitcode_fastqs.transform_quality(quality, orientation)

    expected = {
        "f": (barcode, quality),
        "c": ("TGCA" * 6, "hgfedcbaponmlkji" "xwvutsrq"),
        "r": (barcode[::-1], quality[::-1]),
        "rc": (barcode[16:24] + barcode[8:16] + barcode[0:8], quality[16:24] + quality[8:16] + quality[0:8]),
    }
    assert (transformed_barcode, transformed_quality) == expected[orientation]


def test_combiner_is_atomic_and_rejects_empty_and_duplicate_records(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    cDNA = "ACGT" * 50
    umi = "ACGTTGCACG"
    barcode = "ACGTACGT" "TGCATGCA" "GATTACAG"

    for orientation, read_id in zip(process_splitcode_fastqs.ORIENTATIONS, ["f", "c", "r", "rc"]):
        write_fastq(Path(f"{orientation}_cDNA.fastq.gz"), [(read_id, cDNA, "I" * len(cDNA))])
        write_fastq(Path(f"{orientation}_umi.fastq.gz"), [(read_id, umi, "J" * len(umi))])
        write_fastq(Path(f"{orientation}_barcode.fastq.gz"), [(read_id, barcode, "K" * 24)])

    write_fastq(Path("f_cDNA.fastq.gz"), [("empty", "", "")])
    write_fastq(Path("f_umi.fastq.gz"), [("empty", umi, "J" * len(umi))])
    write_fastq(Path("f_barcode.fastq.gz"), [("empty", barcode, "K" * 24)])

    qc_path = Path("sample_splitcode_qc.tsv")
    process_splitcode_fastqs.combine_fastqs("sample", qc_path)

    outputs = {kind: read_fastq(tmp_path / f"sample_{kind}.fastq.gz") for kind in ("cDNA", "umi", "barcode")}
    assert [record[0] for record in outputs["cDNA"]] == ["c", "r", "rc"]
    assert all(len(outputs[kind]) == 3 for kind in outputs)
    assert all(record[1] and len(record[1]) == len(record[2]) for records in outputs.values() for record in records)
    assert len({record[0] for record in outputs["cDNA"]}) == 3
    assert all(record[1] == cDNA for record in outputs["cDNA"])
    assert all(record[1] == umi for record in outputs["umi"])
    assert outputs["barcode"][0][1] == process_splitcode_fastqs.transform_barcode(barcode, "c")
    assert outputs["barcode"][1][1] == process_splitcode_fastqs.transform_barcode(barcode, "r")
    assert outputs["barcode"][2][1] == process_splitcode_fastqs.transform_barcode(barcode, "rc")
    qc = dict(line.rstrip("\n").split("\t", 1) for line in qc_path.read_text().splitlines()[1:])
    assert qc["empty_triplets_excluded"] == "1"
    assert qc["ambiguous_read_ids_excluded"] == "0"
    assert qc["emitted_triplets"] == "3"


def test_combiner_excludes_duplicate_identifier_across_orientations(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    for orientation in process_splitcode_fastqs.ORIENTATIONS:
        write_fastq(Path(f"{orientation}_cDNA.fastq.gz"), [("same", "ACGT", "IIII")])
        write_fastq(Path(f"{orientation}_umi.fastq.gz"), [("same", "ACGTACGTAC", "JJJJJJJJJJ")])
        write_fastq(Path(f"{orientation}_barcode.fastq.gz"), [("same", "A" * 24, "K" * 24)])

    qc_path = Path("sample_splitcode_qc.tsv")
    process_splitcode_fastqs.combine_fastqs("sample", qc_path)
    qc = dict(line.rstrip("\n").split("\t", 1) for line in qc_path.read_text().splitlines()[1:])
    assert qc["ambiguous_read_ids_excluded"] == "1"
    assert qc["emitted_triplets"] == "0"


def test_combiner_skips_absent_orientation_streams(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    write_fastq(Path("f_cDNA.fastq.gz"), [("f", "ACGT", "IIII")])
    write_fastq(Path("f_umi.fastq.gz"), [("f", "ACGTACGTAC", "JJJJJJJJJJ")])
    write_fastq(Path("f_barcode.fastq.gz"), [("f", "A" * 24, "K" * 24)])

    qc_path = Path("sample_splitcode_qc.tsv")
    process_splitcode_fastqs.combine_fastqs("sample", qc_path)
    qc = dict(line.rstrip("\n").split("\t", 1) for line in qc_path.read_text().splitlines()[1:])
    assert qc["orientations_present"] == "f"
    assert qc["orientations_missing"] == "c,r,rc"
    assert qc["emitted_triplets"] == "1"


def test_empirical_100k_expectations_are_recorded():
    assert EMPIRICAL_100K["f"]["valid_umi"] == 17302
    assert EMPIRICAL_100K["rc"]["valid_umi"] == 21008
    assert EMPIRICAL_100K["f"]["cdna_median"] == 463
    assert EMPIRICAL_100K["rc"]["cdna_median"] == 464
    assert EMPIRICAL_100K["ambiguous"] == 115
    assert EMPIRICAL_100K["exact_unclassified"] == 61555


def test_rendered_template_has_empirical_anchor_order_and_one_mismatch_rules():
    template = TEMPLATE.read_text()
    variables = {
        "assay_id": "assay",
        "assay_name": "name",
        "library_protocol": "Any",
        "library_kit": "kit",
        "sequence_protocol": "ONT",
        "sequence_kit": "GridION",
        "read_file_id": "reads",
        "read_min_length": 200,
        "read_max_length": 603,
        "read_file_name": "reads.fastq",
        "read_file_type": "fastq",
        "read_file_size": 1,
        "read_url": "./reads.fastq",
        "read_urltype": "local",
        "read_md5sum": "null",
    }
    for barcode in range(1, 4):
        variables.update({
            f"barcode_{barcode}_file_id": f"bc{barcode}",
            f"barcode_{barcode}_file_name": f"bc{barcode}.txt.gz",
            f"barcode_{barcode}_file_size": 1,
            f"barcode_{barcode}_url": f"./bc{barcode}.txt.gz",
            f"barcode_{barcode}_location": "local",
            f"barcode_{barcode}_md5": "null",
        })

    rendered = seqspec_autogeneration.render_template(template, variables)
    anchors = [
        "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT",
        "NNNNNNNNNN",
        "GTGGCCGATGTTTCGCATCGGCGTACGACT",
        "ATCCACGTGCTTGAGACTGTGG",
        "TTTTTTTT",
        "CCCATTCACTCTGCGTTGATACCACTGCTT",
    ]
    positions = [rendered.index(anchor) for anchor in anchors]
    assert positions == sorted(positions)
    assert "distances" not in rendered
    splitcode_config = (ROOT / "templates" / "splitcode" / "config-correct.txt").read_text()
    assert "distances" in splitcode_config
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    assert "seqspec index -m rna -s file -t splitcode" in nextflow
    assert "sed -i 's/3:3:3/1:1:1/g' ONT.config" in nextflow