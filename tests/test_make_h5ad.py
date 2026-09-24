import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "make_h5ad.py"
MODULE_SPEC = importlib.util.spec_from_file_location("make_h5ad", SCRIPT)
make_h5ad = importlib.util.module_from_spec(MODULE_SPEC)
MODULE_SPEC.loader.exec_module(make_h5ad)

anndata = pytest.importorskip("anndata")
np = pytest.importorskip("numpy")
mmio = pytest.importorskip("scipy.io")


def write_inputs(tmp_path):
    matrix = tmp_path / "count.mtx"
    barcodes = tmp_path / "count.barcodes.txt"
    features = tmp_path / "count.genes.txt"
    t2g = tmp_path / "reference.t2g"
    mmio.mmwrite(matrix, np.array([[1, 0], [0, 2], [3, 0]], dtype=np.int64))
    barcodes.write_text("ACGTACGTGGGGCCCCTTTTAAAA\nTGCATGCACCCCGGGGAAAATTTT\n")
    features.write_text("tx-1\ntx-2\ntx-3\n")
    t2g.write_text("tx-1\tgene-1\ntx-2\tgene-1\n")
    return matrix, barcodes, features, t2g


def test_transcript_h5ad_is_cell_by_transcript_with_metadata(tmp_path):
    matrix, barcodes, features, t2g = write_inputs(tmp_path)
    output = tmp_path / "sample.transcript.h5ad"

    make_h5ad.convert_matrix(
        matrix_path=matrix,
        barcode_path=barcodes,
        feature_path=features,
        output_path=output,
        feature_type="transcript",
        sample="sample",
        genome="mm39",
        read_type="CDNA",
        entity="nucleus",
        t2g_path=t2g,
    )

    adata = anndata.read_h5ad(output)
    assert adata.shape == (2, 3)
    assert list(adata.obs_names) == ["ACGTACGTGGGGCCCCTTTTAAAA", "TGCATGCACCCCGGGGAAAATTTT"]
    assert list(adata.obs["barcode_1"]) == ["ACGTACGT", "TGCATGCA"]
    assert list(adata.var_names) == ["tx-1", "tx-2", "tx-3"]
    np.testing.assert_array_equal(adata.X.toarray(), [[1, 0, 3], [0, 2, 0]])
    assert list(adata.var["gene_id"]) == ["gene-1", "gene-1", ""]
    assert list(adata.obs["total_counts"]) == [4, 2]
    assert list(adata.obs["n_features"]) == [2, 1]
    assert adata.uns["dogme"]["entity"] == "nucleus"
    assert adata.uns["dogme"]["feature_type"] == "transcript"
    assert adata.uns["dogme"]["unmapped_feature_count"] == 1


def test_matrix_dimensions_must_match_identifiers(tmp_path):
    matrix, barcodes, features, t2g = write_inputs(tmp_path)
    features.write_text("tx-1\ntx-2\n")

    with pytest.raises(ValueError, match="does not match"):
        make_h5ad.convert_matrix(
            matrix_path=matrix,
            barcode_path=barcodes,
            feature_path=features,
            output_path=tmp_path / "invalid.h5ad",
            feature_type="transcript",
            sample="sample",
            genome="mm39",
            read_type="CDNA",
            entity="cell",
            t2g_path=t2g,
        )


def test_barcode_by_feature_matrix_is_not_transposed_twice(tmp_path):
    matrix, barcodes, features, t2g = write_inputs(tmp_path)
    mmio.mmwrite(matrix, np.array([[1, 0, 3], [0, 2, 0]], dtype=np.int64))
    output = tmp_path / "sample.transcript.h5ad"

    make_h5ad.convert_matrix(
        matrix_path=matrix,
        barcode_path=barcodes,
        feature_path=features,
        output_path=output,
        feature_type="transcript",
        sample="sample",
        genome="mm39",
        read_type="CDNA",
        entity="cell",
        t2g_path=t2g,
    )

    adata = anndata.read_h5ad(output)
    assert adata.shape == (2, 3)
    np.testing.assert_array_equal(adata.X.toarray(), [[1, 0, 3], [0, 2, 0]])


def test_qc_reports_h5ad_shapes_and_barcode_umi_thresholds(tmp_path):
    matrix, barcodes, features, t2g = write_inputs(tmp_path)
    mmio.mmwrite(matrix, np.array([[101, 0], [0, 2], [3, 0]], dtype=np.int64))
    gene_h5ad = tmp_path / "sample.gene.h5ad"
    transcript_h5ad = tmp_path / "sample.transcript.h5ad"
    qc_output = tmp_path / "sample.single_cell_qc.tsv"

    make_h5ad.convert_matrix(
        matrix, barcodes, features, gene_h5ad, "gene", "sample", "mm39", "CDNA", "cell"
    )
    make_h5ad.convert_matrix(
        matrix, barcodes, features, transcript_h5ad, "transcript", "sample", "mm39", "CDNA", "cell", t2g
    )
    make_h5ad.write_qc(gene_h5ad, transcript_h5ad, qc_output)

    qc = dict(line.split("\t", 1) for line in qc_output.read_text().splitlines()[1:])
    assert qc["gene_matrix_barcodes"] == "2"
    assert qc["gene_matrix_features"] == "3"
    assert qc["transcript_matrix_barcodes"] == "2"
    assert qc["transcript_matrix_features"] == "3"
    assert qc["barcodes_total"] == "2"
    assert qc["barcodes_with_gt_100_umi"] == "1"
    assert qc["barcodes_with_gt_200_umi"] == "0"


def test_qc_cli_does_not_require_matrix_conversion_arguments():
    args = make_h5ad._build_parser().parse_args([
        "--gene-h5ad", "sample.gene.h5ad",
        "--transcript-h5ad", "sample.transcript.h5ad",
        "--qc-output", "sample.single_cell_qc.tsv",
    ])

    assert args.gene_h5ad == Path("sample.gene.h5ad")
    assert args.transcript_h5ad == Path("sample.transcript.h5ad")
