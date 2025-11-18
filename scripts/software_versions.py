#!/usr/bin/env python3
"""
software_versions.py

Collects versions of the main tools used in the pipeline and writes them to a
softwareVersion.txt file, similar to the old shell-based softwareVTask.
"""

import argparse
import pathlib
import subprocess
from typing import List


def run_command(cmd: List[str]) -> str:
    """
    Run a command and return combined stdout/stderr as a single stripped string.

    On failure, return a descriptive ERROR/COMMAND_NOT_FOUND line instead of
    raising, so the pipeline can still complete and the file is still written.
    """
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            check=True,
        )
        return result.stdout.strip()
    except FileNotFoundError:
        return f"COMMAND_NOT_FOUND: {' '.join(cmd)}"
    except subprocess.CalledProcessError as e:
        out = (e.stdout or "").strip()
        if not out:
            out = str(e)
        return f"ERROR ({' '.join(cmd)}): {out}"


def get_samtools_line(full_output: str) -> str:
    """
    Mimic: samtools version | grep samtools

    Try to return the line containing 'samtools'; fall back to the first line
    or the full output if parsing fails.
    """
    if full_output.startswith("ERROR") or full_output.startswith("COMMAND_NOT_FOUND"):
        return full_output

    lines = [l for l in full_output.splitlines() if l.strip()]
    if not lines:
        return "samtools (no output)"

    for line in lines:
        if "samtools" in line.lower():
            return line

    # Fallback: just return the first non-empty line
    return lines[0]


def list_dorado_models(model_path: str) -> list[str]:
    """
    List entries under the dorado model directory, like the old:
      for folder in "${modelPath}"/*; do basename "$folder"
    """
    p = pathlib.Path(model_path)
    if not p.exists():
        return [f"MODEL_PATH_NOT_FOUND: {model_path}"]

    children = []
    for child in sorted(p.iterdir()):
        children.append(child.name)
    return children


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Write software version information to a file."
    )
    parser.add_argument(
        "--version",
        required=True,
        help="Dogme pipeline version string (used as 'dogme <version>').",
    )
    parser.add_argument(
        "--read-type",
        required=True,
        help="Read type: DNA, RNA, or CDNA (controls which tools are recorded).",
    )
    parser.add_argument(
        "--model-path",
        required=True,
        help="Directory where dorado models are stored.",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output path for the softwareVersion text file.",
    )
    parser.add_argument(
        "--sample",
        required=False,
        help="Sample name (optional, for bookkeeping; not used in file content).",
    )

    args = parser.parse_args()

    read_type = args.read_type.upper()
    out_path = pathlib.Path(args.output)

    lines: list[str] = []

    # 1) dogme version (first line)
    lines.append(f"dogme {args.version}")

    # 2) dorado version (original: dorado -v 2>&1, then echo 'dorado $doradoV')
    dorado_v = run_command(["dorado", "-v"])
    lines.append(f"dorado {dorado_v}")

    # 3) samtools version (original used: samtools version | grep samtools)
    samtools_raw = run_command(["samtools", "version"])
    lines.append(get_samtools_line(samtools_raw))

    # 4) minimap2 version (original: minimap2 --version 2>&1, echo 'minimap2 $minimap2V')
    minimap2_v = run_command(["minimap2", "--version"])
    lines.append(f"minimap2 {minimap2_v}")

    # 5) modkit if DNA or RNA
    if read_type in ("DNA", "RNA"):
        modkit_v = run_command(["modkit", "--version"])
        lines.append(modkit_v)

    # 6) kallisto + bustools if CDNA or RNA
    if read_type in ("CDNA", "RNA"):
        kallisto_v = run_command(["kallisto", "version"])
        lines.append(kallisto_v)

        bustools_v = run_command(["bustools", "version"])
        lines.append(bustools_v)

    # 7) Dorado models used
    lines.append("Dorado Models Used: ")
    for model_name in list_dorado_models(args.model_path):
        lines.append(model_name)

    # Write file
    out_path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()

    