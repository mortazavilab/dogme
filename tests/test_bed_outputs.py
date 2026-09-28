import shutil
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FILTER_SCRIPT = ROOT / "scripts" / "filterbed.py"


def test_filterbed_writes_plain_bed_for_downstream_bgzip(tmp_path):
    bgzip = shutil.which("bgzip")
    if bgzip is None:
        pytest.skip("bgzip is unavailable")

    input_bed = tmp_path / "sample.bed.gz"
    output_bed = tmp_path / "sample.filtered-3-5.bed"
    rows = [
        "chr2\t20\t30\t.\t0\t+\t.\t.\t.\t.\t6\t10\n",
        "chr1\t2\t4\t.\t0\t+\t.\t.\t.\t.\t4\t10\n",
    ]

    input_plain = tmp_path / "sample.bed"
    with input_plain.open("wt") as handle:
        handle.writelines(rows)
    with input_bed.open("wb") as compressed:
        subprocess.run([bgzip, "-c", str(input_plain)], check=True, stdout=compressed)

    subprocess.run(
        [sys.executable, str(FILTER_SCRIPT), "3", "5", str(input_bed), str(output_bed)],
        check=True,
        capture_output=True,
        text=True,
    )

    expected_row = rows[0].split("\t")
    expected_row[4] = "6"
    assert output_bed.read_text() == "\t".join(expected_row)
    assert output_bed.read_bytes()[:2] != b"\x1f\x8b"


def test_bed_processes_declare_bgzip_and_tabix_outputs():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()

    assert nextflow.count("bgzip -@ ${task.cpus}") >= 4
    assert nextflow.count("tabix -f -p bed") >= 4
    assert 'path "*.bed.gz.tbi", emit: indexes' in nextflow
    assert 'path "${params.sample}.${genomeName}.m6Aopen.bed.gz", emit: beds' in nextflow
    assert "params.readType == 'CDNA' && singleCellEnabled" in nextflow


def test_fastq_processes_use_task_allocated_bgzip_threads():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    extract_fastq = nextflow.split("process extractfastqTask", 1)[1].split("process demuxExtractfastqTask", 1)[0]
    demux_extract_fastq = nextflow.split("process demuxExtractfastqTask", 1)[1].split("process demuxGenerateSeqspecTask", 1)[0]
    splitcode = nextflow.split("process splitcodeTask", 1)[1].split("process singleCellKallistoTask", 1)[0]

    assert "samtools fastq --threads ${task.cpus}" in extract_fastq
    assert "bgzip -@ ${task.cpus} -f ${params.sample}.fastq" in extract_fastq
    assert "samtools fastq --threads ${task.cpus}" in demux_extract_fastq
    assert "bgzip -@ ${task.cpus} -f ${sampleName}.fastq" in demux_extract_fastq
    assert 'bgzip -@ ${task.cpus} -f "${params.sample}_barcode.fastq"' in splitcode