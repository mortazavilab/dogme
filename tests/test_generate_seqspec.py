import importlib.util
import gzip
import json
import shutil
import subprocess
import sys
from pathlib import Path


import pytest

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = ROOT / "tests" / "fixtures"
SCRIPT = ROOT / "scripts" / "generate_seqspec.py"
AUTOGENERATION = ROOT / "scripts" / "seqspec_autogeneration.py"
module_spec = importlib.util.spec_from_file_location("seqspec_autogeneration", AUTOGENERATION)
seqspec_autogeneration = importlib.util.module_from_spec(module_spec)
module_spec.loader.exec_module(seqspec_autogeneration)
load_overrides = seqspec_autogeneration.load_overrides
load_single_cell_kit_variables = seqspec_autogeneration.load_single_cell_kit_variables
measure_fastq = seqspec_autogeneration.measure_fastq
select_template = seqspec_autogeneration.select_template
build_variables = seqspec_autogeneration.build_variables
render_template = seqspec_autogeneration.render_template


@pytest.mark.parametrize(
    ("read_type", "single_cell", "expected"),
    [
        ("RNA", False, "ont-bulk-drna.yaml.j2"),
        ("CDNA", False, "ont-bulk-cdna.yaml.j2"),
        ("DNA", False, "ont-bulk-gdna.yaml.j2"),
    ],
)
def test_select_builtin_bulk_template(read_type, single_cell, expected):
    assert select_template(read_type, single_cell, ROOT / "templates" / "seqspec").name == expected


def test_select_builtin_single_cell_template():
    assert select_template("CDNA", True, ROOT / "templates" / "seqspec").name == "parse-evercode-wt-mega-v2-nanopore.yaml.j2"


def test_unsupported_template_combination_names_inputs():
    with pytest.raises(ValueError, match="readType=DNA singleCell=true"):
        select_template("DNA", True, ROOT / "templates" / "seqspec")


@pytest.mark.parametrize("kit", ["parse-wt-v2", "parse-wt-mega-v2"])
def test_parse_kit_defaults_include_onlist_metadata(kit):
    variables = load_single_cell_kit_variables(ROOT / "templates" / "seqspec", kit)

    assert variables["assay_id"]
    assert variables["library_kit"]
    for barcode in range(1, 4):
        assert variables[f"barcode_{barcode}_file_id"]
        assert variables[f"barcode_{barcode}_file_name"].endswith(".txt.gz")
        assert variables[f"barcode_{barcode}_url"].startswith("https://")
        assert variables[f"barcode_{barcode}_md5"]


def test_parse_kit_defaults_reject_unknown_kit():
    with pytest.raises(ValueError, match="unsupported singleCellKit=unknown"):
        load_single_cell_kit_variables(ROOT / "templates" / "seqspec", "unknown")


@pytest.mark.parametrize("kit", ["parse-wt-v2", "parse-wt-mega-v2"])
def test_parse_kit_defaults_render_every_template_placeholder(tmp_path, kit):
    fastq = tmp_path / "sample.fastq"
    fastq.write_text("@read\nACGT\n+\n!!!!\n")
    template = (ROOT / "templates" / "parse-evercode-wt-mega-v2-nanopore.yaml.j2").read_text()
    variables = build_variables(
        fastq,
        load_single_cell_kit_variables(ROOT / "templates" / "seqspec", kit),
    )

    rendered = render_template(template, variables)

    assert "{{" not in rendered
    assert f"assay_id: {variables['assay_id']}" in rendered


def test_measure_fastq_and_merge_overrides(tmp_path):
    fastq = tmp_path / "sample.fastq"
    fastq.write_text("@one\nACGT\n+\n!!!!\n@two\nAC\n+\n!!\n")
    measured = measure_fastq(fastq)
    assert measured["read_min_length"] == 2
    assert measured["read_max_length"] == 4
    assert measured["read_file_size"] == fastq.stat().st_size
    assert measured["read_md5sum"]

    variables = tmp_path / "variables.json"
    variables.write_text('{"read_max_length": 99, "library_kit": "override"}')
    overrides = load_overrides(variables)
    assert overrides == {"read_max_length": 99, "library_kit": "override"}
    merged = build_variables(fastq, overrides)
    assert merged["read_max_length"] == 99
    assert merged["read_min_length"] == 2
    assert merged["library_kit"] == "override"


def test_measure_gzipped_fastq(tmp_path):
    fastq = tmp_path / "sample.fastq.gz"
    with gzip.open(fastq, "wt") as handle:
        handle.write("@one\nACGT\n+\n!!!!\n")

    measured = measure_fastq(fastq)

    assert measured["read_file_name"] == "sample.fastq.gz"
    assert measured["read_min_length"] == 4
    assert measured["read_max_length"] == 4
    assert measured["read_file_size"] == fastq.stat().st_size


def test_render_template_replaces_values_and_arithmetic():
    assert render_template(
        "name={{ name }} max={{ read_max_length - 154 }} empty={{ missing_value }}",
        {"name": "sample", "read_max_length": 200, "missing_value": None},
    ) == "name=sample max=46 empty=null"


def test_render_template_reports_missing_variable():
    with pytest.raises(KeyError, match="missing seqspec template variable: name"):
        render_template("{{ name }}", {})


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


def test_seqspec_031_version_is_accepted_despite_nonzero_exit(monkeypatch):
    module_spec = importlib.util.spec_from_file_location("generate_seqspec", SCRIPT)
    generate_seqspec = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(generate_seqspec)

    monkeypatch.setattr(
        generate_seqspec.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(args[0], 1, "seqspec 0.3.1\n"),
    )

    assert generate_seqspec.get_seqspec_version() == (0, 3, 1)


@pytest.mark.parametrize(
    ("version", "expected_calls"),
    [
        ((0, 3, 1), [("format", "rendered.yaml", "-o", "output.yaml"), ("check", "output.yaml")]),
        (
            (0, 4, 0),
            [
                ("upgrade", "rendered.yaml", "-o", "output.yaml"),
                ("format", "output.yaml", "-o", "output.yaml"),
                ("check", "output.yaml"),
            ],
        ),
    ],
)
def test_finalize_seqspec_uses_version_appropriate_commands(monkeypatch, version, expected_calls):
    module_spec = importlib.util.spec_from_file_location("generate_seqspec", SCRIPT)
    generate_seqspec = importlib.util.module_from_spec(module_spec)
    module_spec.loader.exec_module(generate_seqspec)
    calls = []

    monkeypatch.setattr(generate_seqspec, "run_seqspec", lambda *args: calls.append(args))

    generate_seqspec.finalize_seqspec(Path("rendered.yaml"), Path("output.yaml"), version)

    assert calls == expected_calls


def test_generate_seqspec_with_image_runtime(tmp_path):
    if shutil.which("seqspec") is None:
        pytest.skip("seqspec is unavailable; run this test inside the DOGME image")

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