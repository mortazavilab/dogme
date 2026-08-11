#!/usr/bin/env python3
"""Render and validate a seqspec artifact for a generated FASTQ."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path
import re

try:
    from seqspec_autogeneration import (
        build_variables,
        load_overrides,
        load_single_cell_kit_variables,
        render_template,
        select_template,
    )
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    from seqspec_autogeneration import (
        build_variables,
        load_overrides,
        load_single_cell_kit_variables,
        render_template,
        select_template,
    )

SEQSPEC_VERSION_PATTERN = re.compile(r"seqspec\s+(\d+)\.(\d+)\.(\d+)", re.IGNORECASE)


def run_seqspec(*args: str) -> str:
    command = ["seqspec", *args]
    try:
        result = subprocess.run(
            command,
            check=True,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "seqspec is not available in the Docker/Apptainer image. "
            "Use a DOGME image containing seqspec 0.4.0."
        ) from exc
    except subprocess.CalledProcessError as exc:
        output = (exc.stdout or "").strip()
        raise RuntimeError(f"seqspec {' '.join(args)} failed:\n{output}") from exc
    return result.stdout


def get_seqspec_version() -> tuple[int, int, int]:
    try:
        result = subprocess.run(
            ["seqspec", "--version"],
            check=False,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
    except FileNotFoundError as exc:
        raise RuntimeError(
            "seqspec is not available in the Docker/Apptainer image. "
            "Install seqspec 0.3.1 or newer."
        ) from exc

    output = result.stdout or ""
    match = SEQSPEC_VERSION_PATTERN.search(output)
    if match is None:
        raise RuntimeError(f"Could not determine seqspec version:\n{output.strip()}")
    return tuple(int(value) for value in match.groups())


def finalize_seqspec(rendered_path: Path, output_path: Path, version: tuple[int, int, int]) -> None:
    if version >= (0, 4, 0):
        run_seqspec("upgrade", str(rendered_path), "-o", str(output_path))
        run_seqspec("format", str(output_path), "-o", str(output_path))
    elif version >= (0, 3, 0):
        run_seqspec("format", str(rendered_path), "-o", str(output_path))
    else:
        raise RuntimeError(
            f"seqspec {'.'.join(str(value) for value in version)} is unsupported; "
            "install seqspec 0.3.1 or newer."
        )
    run_seqspec("check", str(output_path))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", type=Path)
    parser.add_argument("--variables", type=Path)
    parser.add_argument("--template-dir", type=Path)
    parser.add_argument("--read-type", default="RNA")
    parser.add_argument("--single-cell", action="store_true")
    parser.add_argument("--single-cell-kit")
    parser.add_argument("--no-md5", action="store_true")
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--fastq", required=True, type=Path)
    args = parser.parse_args()

    uses_builtin_template = args.template is None
    if args.template is None:
        if args.template_dir is None:
            parser.error("--template or --template-dir is required")
        try:
            args.template = select_template(args.read_type, args.single_cell, args.template_dir)
        except (FileNotFoundError, ValueError) as exc:
            parser.error(str(exc))
    if not args.template.is_file():
        parser.error(f"seqspec template does not exist: {args.template}")
    if not args.fastq.is_file():
        parser.error(f"FASTQ does not exist: {args.fastq}")

    if args.single_cell_kit and not args.single_cell:
        parser.error("--single-cell-kit requires --single-cell")
    try:
        selected_kit = args.single_cell_kit
        if uses_builtin_template and args.single_cell and selected_kit is None:
            selected_kit = "parse-wt-mega-v2"
        kit_variables = load_single_cell_kit_variables(args.template_dir, selected_kit)
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))
    variables = build_variables(args.fastq, kit_variables, load_overrides(args.variables))
    if args.no_md5:
        variables["read_md5sum"] = None


    variables.setdefault("read_file_name", args.fastq.name)
    try:
        rendered = render_template(args.template.read_text(), variables)
    except (KeyError, TypeError, ValueError) as exc:
        raise RuntimeError(f"Could not render seqspec template: {exc}") from exc
    rendered_path = args.output.with_suffix(".rendered.yaml")
    rendered_path.write_text(rendered)

    finalize_seqspec(rendered_path, args.output, get_seqspec_version())
    args.output.with_suffix(".variables.json").write_text(
        json.dumps(variables, indent=2, sort_keys=True) + "\n"
    )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)