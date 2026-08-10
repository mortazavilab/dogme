#!/usr/bin/env python3
"""Render and validate a seqspec artifact for a generated FASTQ."""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

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


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--template", required=True, type=Path)
    parser.add_argument("--variables", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--fastq", required=True, type=Path)
    args = parser.parse_args()

    if not args.template.is_file():
        parser.error(f"seqspec template does not exist: {args.template}")
    if not args.fastq.is_file():
        parser.error(f"FASTQ does not exist: {args.fastq}")

    try:
        variables = json.loads(args.variables.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Could not read seqspec render variables: {exc}") from exc

    try:
        from jinja2 import Template
    except ModuleNotFoundError as exc:
        raise RuntimeError(
            "Jinja2 is not available in the DOGME Docker/Apptainer image; "
            "the image must provide Jinja2 for seqspec template rendering."
        ) from exc

    variables.setdefault("read_file_name", args.fastq.name)
    rendered = Template(args.template.read_text()).render(**variables)
    rendered_path = args.output.with_suffix(".rendered.yaml")
    rendered_path.write_text(rendered)

    # Verify the image exposes the requested CLI before invoking it. This also
    # gives a useful failure when a host runs the workflow without the image.
    run_seqspec("--version")
    run_seqspec("check", "--help")
    run_seqspec("upgrade", "--help")
    run_seqspec("format", "--help")

    # seqspec 0.4.0 writes upgraded/ formatted specs to the requested path.
    run_seqspec("upgrade", str(rendered_path), "-o", str(args.output))
    run_seqspec("format", str(args.output), "-o", str(args.output))
    run_seqspec("check", str(args.output))
    shutil.copyfile(args.variables, args.output.with_suffix(".variables.json"))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except RuntimeError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(1)