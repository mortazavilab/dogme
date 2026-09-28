#!/usr/bin/env python3
"""Convert a bustools Matrix Market count matrix to an AnnData H5AD file."""

from __future__ import annotations

import argparse
from pathlib import Path


def _load_dependencies():
    try:
        import anndata
        import h5py
        import numpy as np
        import pandas as pd
        from scipy import io, sparse
    except ImportError as exc:
        raise RuntimeError(
            "H5AD generation requires anndata, h5py, numpy, pandas, and scipy "
            "in the DOGME execution image. Install requirements-h5ad.txt."
        ) from exc
    return anndata, io, np, pd, sparse


def _read_identifiers(path: Path) -> list[str]:
    identifiers = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            identifier = fields[0].strip()
            if not identifier:
                raise ValueError(f"Empty feature identifier in {path} at line {line_number}")
            identifiers.append(identifier)
    if not identifiers:
        raise ValueError(f"No identifiers found in {path}")
    if len(set(identifiers)) != len(identifiers):
        raise ValueError(f"Duplicate identifiers found in {path}")
    return identifiers


def _read_t2g(path: Path) -> dict[str, str]:
    mapping = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 2 or not fields[0] or not fields[1]:
                raise ValueError(f"Invalid transcript-to-gene row in {path} at line {line_number}")
            transcript_id, gene_id = fields[0], fields[1]
            if transcript_id in mapping and mapping[transcript_id] != gene_id:
                raise ValueError(f"Conflicting transcript mapping for {transcript_id} in {path}")
            mapping[transcript_id] = gene_id
    return mapping


def _as_integer_matrix(matrix, np, sparse):
    matrix = sparse.csr_matrix(matrix)
    if matrix.data.size and (matrix.data < 0).any():
        raise ValueError("Count matrix contains negative values")
    if matrix.data.size and not np.all(np.equal(matrix.data, np.floor(matrix.data))):
        raise ValueError("Count matrix contains non-integer values")
    return matrix.astype(np.int64, copy=False)


def convert_matrix(
    matrix_path: Path,
    barcode_path: Path,
    feature_path: Path,
    output_path: Path,
    feature_type: str,
    sample: str,
    genome: str,
    read_type: str,
    entity: str,
    t2g_path: Path | None = None,
) -> None:
    """Convert a bustools count matrix into cell-by-feature AnnData."""
    anndata, mmio, np, pd, sparse = _load_dependencies()
    for path in (matrix_path, barcode_path, feature_path):
        if not path.exists():
            raise FileNotFoundError(f"Required count input does not exist: {path}")
    if feature_type not in {"gene", "transcript"}:
        raise ValueError(f"Unsupported feature type: {feature_type}")
    if entity not in {"cell", "nucleus"}:
        raise ValueError("entity must be either 'cell' or 'nucleus'")

    barcodes = _read_identifiers(barcode_path)
    features = _read_identifiers(feature_path)
    source_counts = _as_integer_matrix(mmio.mmread(matrix_path), np, sparse)
    feature_by_barcode_shape = (len(features), len(barcodes))
    barcode_by_feature_shape = (len(barcodes), len(features))
    if source_counts.shape == feature_by_barcode_shape:
        counts = source_counts.transpose().tocsr()
    elif source_counts.shape == barcode_by_feature_shape:
        counts = source_counts
    else:
        raise ValueError(
            f"Matrix shape {source_counts.shape} does not match either "
            f"{len(features)} features x {len(barcodes)} barcodes or "
            f"{len(barcodes)} barcodes x {len(features)} features"
        )
    obs = pd.DataFrame(index=pd.Index(barcodes, name="barcode"))
    obs["barcode_1"] = [barcode[:8] for barcode in barcodes]
    var = pd.DataFrame(index=pd.Index(features, name=feature_type))
    obs["total_counts"] = np.asarray(counts.sum(axis=1)).ravel().astype(np.int64)
    obs["n_features"] = np.asarray((counts > 0).sum(axis=1)).ravel().astype(np.int64)
    var["total_counts"] = np.asarray(counts.sum(axis=0)).ravel().astype(np.int64)
    var["n_cells"] = np.asarray((counts > 0).sum(axis=0)).ravel().astype(np.int64)

    metadata = {
        "sample": sample,
        "genome": genome,
        "read_type": read_type,
        "entity": entity,
        "feature_type": feature_type,
        "matrix_orientation": "cells_by_features",
        "count_type": "raw_integer_counts",
        "source_matrix": str(matrix_path.name),
        "source_barcodes": str(barcode_path.name),
        "source_features": str(feature_path.name),
    }
    if feature_type == "transcript" and t2g_path is not None:
        if not t2g_path.exists():
            raise FileNotFoundError(f"Transcript-to-gene mapping does not exist: {t2g_path}")
        transcript_to_gene = _read_t2g(t2g_path)
        var["gene_id"] = [transcript_to_gene.get(transcript_id, "") for transcript_id in features]
        metadata["source_t2g"] = str(t2g_path.name)
        metadata["unmapped_feature_count"] = sum(
            transcript_id not in transcript_to_gene for transcript_id in features
        )

    adata = anndata.AnnData(X=counts, obs=obs, var=var)
    adata.uns["dogme"] = metadata
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path, compression="gzip")


def write_qc(gene_h5ad_path: Path, transcript_h5ad_path: Path, output_path: Path) -> None:
    """Write matrix dimensions and gene-level UMI barcode thresholds."""
    anndata, _, np, _, _ = _load_dependencies()
    for path in (gene_h5ad_path, transcript_h5ad_path):
        if not path.exists():
            raise FileNotFoundError(f"Required H5AD input does not exist: {path}")

    gene_data = anndata.read_h5ad(gene_h5ad_path, backed="r")
    transcript_data = anndata.read_h5ad(transcript_h5ad_path, backed="r")
    try:
        if "total_counts" not in gene_data.obs:
            raise ValueError(f"Gene H5AD is missing obs['total_counts']: {gene_h5ad_path}")

        umi_counts = np.asarray(gene_data.obs["total_counts"], dtype=np.int64)
        metrics: list[tuple[str, int | str]] = [
            ("gene_matrix_barcodes", gene_data.n_obs),
            ("gene_matrix_features", gene_data.n_vars),
            ("transcript_matrix_barcodes", transcript_data.n_obs),
            ("transcript_matrix_features", transcript_data.n_vars),
            ("barcode_umi_source", "gene_h5ad_obs_total_counts"),
            ("barcodes_total", len(umi_counts)),
        ]
        for threshold in (100, 200, 300, 400, 500, 1000, 2000, 3000, 4000, 10000, 20000):
            metrics.append((f"barcodes_with_gt_{threshold}_umi", int((umi_counts > threshold).sum())))
    finally:
        gene_data.file.close()
        transcript_data.file.close()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as output:
        output.write("metric\tvalue\n")
        for metric, value in metrics:
            output.write(f"{metric}\t{value}\n")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path)
    parser.add_argument("--barcodes", type=Path)
    parser.add_argument("--features", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--feature-type", choices=("gene", "transcript"))
    parser.add_argument("--sample")
    parser.add_argument("--genome")
    parser.add_argument("--read-type")
    parser.add_argument("--entity", choices=("cell", "nucleus"), default="cell")
    parser.add_argument("--t2g", type=Path)
    parser.add_argument("--gene-h5ad", type=Path)
    parser.add_argument("--transcript-h5ad", type=Path)
    parser.add_argument("--qc-output", type=Path)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    if args.qc_output is not None:
        if args.gene_h5ad is None or args.transcript_h5ad is None:
            raise ValueError("--qc-output requires --gene-h5ad and --transcript-h5ad")
        write_qc(args.gene_h5ad, args.transcript_h5ad, args.qc_output)
        return
    if args.gene_h5ad is not None or args.transcript_h5ad is not None:
        raise ValueError("--gene-h5ad and --transcript-h5ad require --qc-output")
    required_arguments = {
        "--matrix": args.matrix,
        "--barcodes": args.barcodes,
        "--features": args.features,
        "--output": args.output,
        "--feature-type": args.feature_type,
        "--sample": args.sample,
        "--genome": args.genome,
        "--read-type": args.read_type,
    }
    missing = [name for name, value in required_arguments.items() if value is None]
    if missing:
        raise ValueError(f"H5AD conversion requires: {', '.join(missing)}")
    convert_matrix(
        matrix_path=args.matrix,
        barcode_path=args.barcodes,
        feature_path=args.features,
        output_path=args.output,
        feature_type=args.feature_type,
        sample=args.sample,
        genome=args.genome,
        read_type=args.read_type,
        entity=args.entity,
        t2g_path=args.t2g,
    )


if __name__ == "__main__":
    main()
