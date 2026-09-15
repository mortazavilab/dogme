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
    barcodes.write_text("cell-a\ncell-b\n")
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
    assert list(adata.obs_names) == ["cell-a", "cell-b"]
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
