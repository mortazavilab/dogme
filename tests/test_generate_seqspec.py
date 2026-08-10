import importlib.util
import json
import shutil
import subprocess
import sys
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"
SCRIPT = ROOT / "scripts" / "generate_seqspec.py"


def run_script(*args):
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )


def test_argument_help_and_required_arguments():
    help_result = run_script("--help")
    assert help_result.returncode == 0
    assert "--template" in help_result.stdout

    missing_argument_result = run_script()
    assert missing_argument_result.returncode == 2
    assert "the following arguments are required" in missing_argument_result.stderr


def test_missing_template_error(tmp_path):
    result = run_script(
        "--template",
        str(tmp_path / "missing.yaml.j2"),
        "--variables",
        str(FIXTURES / "synthetic-seqspec-variables.json"),
        "--fastq",
        str(tmp_path / "missing.fastq"),
        "--output",
        str(tmp_path / "sample.seqspec.yaml"),
    )
    assert result.returncode == 2
    assert "seqspec template does not exist" in result.stderr


def test_missing_fastq_error(tmp_path):
    template = tmp_path / "template.yaml.j2"
    template.write_text("placeholder")
    result = run_script(
        "--template",
        str(template),
        "--variables",
        str(tmp_path / "missing.json"),
        "--fastq",
        str(tmp_path / "missing.fastq"),
        "--output",
        str(tmp_path / "sample.seqspec.yaml"),
    )
    assert result.returncode == 2
    assert "FASTQ does not exist" in result.stderr


def test_missing_variables_json_after_fastq_validation(tmp_path):
    template = tmp_path / "template.yaml.j2"
    template.write_text("placeholder")
    fastq = tmp_path / "input.fastq"
    fastq.write_text("@read\nACGT\n+\n!!!!\n")
    result = run_script(
        "--template",
        str(template),
        "--variables",
        str(tmp_path / "missing.json"),
        "--fastq",
        str(fastq),
        "--output",
        str(tmp_path / "sample.seqspec.yaml"),
    )
    assert result.returncode == 1
    assert "Could not read seqspec render variables" in result.stderr


def test_seqspec_not_found_error(monkeypatch):
    module_spec = importlib.util.spec_from_file_location("generate_seqspec", SCRIPT)
    generate_seqspec = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(generate_seqspec)

    def missing_executable(*args, **kwargs):
        raise FileNotFoundError("seqspec")

    monkeypatch.setattr(generate_seqspec.subprocess, "run", missing_executable)
    with pytest.raises(RuntimeError, match="seqspec is not available"):
        generate_seqspec.run_seqspec("--version")


def test_generate_seqspec_with_image_runtime(tmp_path):
    if shutil.which("seqspec") is None:
        pytest.skip("seqspec is unavailable; run this test inside the DOGME image")
    if importlib.util.find_spec("jinja2") is None:
        pytest.skip("Jinja2 is unavailable; run this test inside the DOGME image")

    fastq = tmp_path / "synthetic.fastq"
    fastq.write_text("@read-1\nACGT\n+\n!!!!\n")
    variables = json.loads((FIXTURES / "synthetic-seqspec-variables.json").read_text())
    variables_path = tmp_path / "variables.json"
    variables_path.write_text(json.dumps(variables))
    output = tmp_path / "sample.seqspec.yaml"

    subprocess.run(
        [
            "python3",
            str(ROOT / "scripts" / "generate_seqspec.py"),
            "--template",
            str(FIXTURES / "synthetic-seqspec.yaml.j2"),
            "--variables",
            str(variables_path),
            "--fastq",
            str(fastq),
            "--output",
            str(output),
        ],
        check=True,
        cwd=ROOT,
    )

    assert output.is_file()
    assert output.with_name("sample.seqspec.rendered.yaml").is_file()
    assert output.with_name("sample.seqspec.variables.json").is_file()