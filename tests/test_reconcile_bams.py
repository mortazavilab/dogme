import importlib.util
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "reconcileBams.py"
pytest.importorskip("pysam")
MODULE_SPEC = importlib.util.spec_from_file_location("reconcile_bams", SCRIPT)
reconcile_bams = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(reconcile_bams)


@pytest.mark.parametrize("entity", ("cell", "nucleus"))
def test_recreates_h5ad_outputs_from_reconciled_annotated_bam(tmp_path, monkeypatch, entity):
    calls = []
    monkeypatch.setattr(
        reconcile_bams.subprocess,
        "run",
        lambda command, check: calls.append((command, check)),
    )

    input_bam = "/data/sample.mm39.annotated.single_cell.bam"
    reconciled_bam = str(tmp_path / "sample.mm39.annotated.single_cell.reconciled.bam")
    gtf = str(tmp_path / "combined.gtf")
    reconcile_bams.recreate_h5ad(
        input_bam,
        reconciled_bam,
        gtf,
        str(tmp_path),
        "mm39",
        entity,
    )

    command, check = calls[0]
    assert command[:3] == [
        sys.executable,
        str(ROOT / "scripts" / "make_h5ad.py"),
        "--bam",
    ]
    assert command[2:6] == ["--bam", reconciled_bam, "--gtf", gtf]
    assert command[command.index("--gene-output") + 1] == str(tmp_path / "sample.mm39.gene.h5ad")
    assert command[command.index("--transcript-output") + 1] == str(tmp_path / "sample.mm39.transcript.h5ad")
    assert command[command.index("--qc-output") + 1] == str(tmp_path / "sample.mm39.single_cell_qc.tsv")
    assert command[command.index("--sample") + 1] == "sample"
    assert command[command.index("--genome") + 1] == "mm39"
    assert command[command.index("--entity") + 1] == entity
    assert check is True


def test_recreation_requires_genome_argument(tmp_path, monkeypatch, capsys):
    monkeypatch.setattr(
        sys,
        "argv",
        [
            str(SCRIPT),
            "--bams", str(tmp_path / "sample.bam"),
            "--annotation", str(tmp_path / "reference.gtf"),
            "--out_prefix", "sample",
            "--outdir", str(tmp_path),
            "--recreate-h5ad",
        ],
    )

    with pytest.raises(SystemExit) as error:
        reconcile_bams.main()

    assert error.value.code == 2
    assert "--genome is required when --recreate-h5ad is enabled" in capsys.readouterr().err
