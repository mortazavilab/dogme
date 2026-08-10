"""Select built-in seqspec templates and derive FASTQ metadata."""

from __future__ import annotations

import hashlib
import json
from datetime import date
from pathlib import Path
from typing import Any


def select_template(read_type: str, single_cell: bool, template_dir: Path) -> Path:
    normalized_type = read_type.upper()
    if single_cell:
        if normalized_type != "CDNA":
            raise ValueError(
                f"no seqspec template for readType={normalized_type} singleCell=true"
            )
        filename = "parse-evercode-wt-mega-v2-nanopore.yaml.j2"
    else:
        filenames = {
            "RNA": "ont-bulk-drna.yaml.j2",
            "CDNA": "ont-bulk-cdna.yaml.j2",
            "DNA": "ont-bulk-gdna.yaml.j2",
        }
        filename = filenames.get(normalized_type)
        if filename is None:
            raise ValueError(
                f"no seqspec template for readType={normalized_type} singleCell=false"
            )
    template = template_dir / filename
    if not template.is_file() and template_dir.name == "seqspec":
        template = template_dir.parent / filename
    if not template.is_file():
        raise FileNotFoundError(f"built-in seqspec template does not exist: {template}")
    return template


def measure_fastq(fastq: Path, calculate_md5: bool = True) -> dict[str, Any]:
    minimum: int | None = None
    maximum: int | None = None
    size = fastq.stat().st_size
    digest = hashlib.md5() if calculate_md5 else None
    with fastq.open("rb") as handle:
        while True:
            record = [handle.readline() for _ in range(4)]
            if not record[0]:
                break
            if len(record) != 4 or not record[2].startswith(b"+"):
                raise ValueError(f"invalid FASTQ record in {fastq}")
            length = len(record[1].rstrip(b"\r\n"))
            minimum = length if minimum is None else min(minimum, length)
            maximum = length if maximum is None else max(maximum, length)
            if digest is not None:
                digest.update(b"".join(record))
    if minimum is None:
        raise ValueError(f"FASTQ contains no reads: {fastq}")
    return {
        "date": date.today().strftime("%-d %B %Y"),
        "read_file_name": fastq.name,
        "read_file_id": fastq.name,
        "read_file_type": "fastq",
        "read_file_size": size,
        "read_url": f"./{fastq.name}",
        "read_urltype": "local",
        "read_md5sum": digest.hexdigest() if digest is not None else None,
        "read_min_length": minimum,
        "read_max_length": maximum,
    }


def load_overrides(path: Path | None) -> dict[str, Any]:
    if path is None:
        return {}
    try:
        value = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError) as exc:
        raise RuntimeError(f"Could not read seqspec render variables: {exc}") from exc
    if not isinstance(value, dict):
        raise ValueError(f"seqspec render variables must be a JSON object: {path}")
    return value


def build_variables(fastq: Path, overrides: dict[str, Any] | None = None) -> dict[str, Any]:
    variables = measure_fastq(fastq)
    variables.update(overrides or {})
    return variables