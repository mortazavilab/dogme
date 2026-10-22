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


def test_annotated_bam_h5ad_includes_novel_tags_and_deduplicates_umis(tmp_path):
    pysam = pytest.importorskip("pysam")
    bam_path = tmp_path / "sample.annotated.single_cell.bam"
    gtf_path = tmp_path / "sample_dogme.gtf"
    gene_output = tmp_path / "sample.mm39.gene.h5ad"
    transcript_output = tmp_path / "sample.mm39.transcript.h5ad"
    qc_output = tmp_path / "sample.mm39.single_cell_qc.tsv"
    gtf_path.write_text(
        'chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id "gene-known"; '
        'transcript_id "tx-known"; gene_name "Known"; gene_type "protein_coding"; '
        'transcript_name "Known transcript"; transcript_type "protein_coding";\n'
        'chr1\ttest\texon\t20\t30\t.\t+\t.\tgene_id "gene-unused"; '
        'transcript_id "tx-unused"; gene_name "Unused"; gene_type "lncRNA"; '
        'transcript_name "Unused transcript"; transcript_type "lncRNA";\n'
    )
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 1000}]}
    records = [
        ("read1", "CELL0001", "UMI0000001", "gene-known", "tx-known"),
        ("read1-duplicate", "CELL0001", "UMI0000001", "gene-known", "tx-known"),
        ("read2", "CELL0001", "UMI0000002", "NOVELG000001", "NOVELT0000000001"),
        ("read3", "CELL0002", "UMI0000001", "gene-known", "tx-known"),
    ]
    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for position, (read_id, barcode, umi, gene_id, transcript_id) in enumerate(records):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = read_id
            read.query_sequence = "A" * 10
            read.flag = 0
            read.reference_id = 0
            read.reference_start = position * 10
            read.mapping_quality = 60
            read.cigar = ((0, 10),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 10)
            read.set_tag("CB", barcode)
            read.set_tag("UB", umi)
            read.set_tag("GX", gene_id)
            read.set_tag("TX", transcript_id)
            bam.write(read)

    make_h5ad.convert_annotated_bam(
        bam_path=bam_path,
        gtf_path=gtf_path,
        gene_output_path=gene_output,
        transcript_output_path=transcript_output,
        sample="sample",
        genome="mm39",
        entity="cell",
        qc_output_path=qc_output,
    )

    gene_data = anndata.read_h5ad(gene_output)
    transcript_data = anndata.read_h5ad(transcript_output)
    assert set(gene_data.var_names) == {"gene-known", "gene-unused", "NOVELG000001"}
    assert set(transcript_data.var_names) == {"tx-known", "tx-unused", "NOVELT0000000001"}
    assert list(gene_data.obs_names) == ["CELL0001", "CELL0002"]
    assert list(gene_data.obs["barcode_1"]) == ["CELL0001", "CELL0002"]
    assert gene_data[0, "gene-known"].X[0, 0] == 1
    assert gene_data[0, "NOVELG000001"].X[0, 0] == 1
    assert gene_data[1, "gene-known"].X[0, 0] == 1
    assert transcript_data[0, "tx-known"].X[0, 0] == 1
    assert transcript_data[0, "NOVELT0000000001"].X[0, 0] == 1
    assert transcript_data[1, "tx-known"].X[0, 0] == 1
    assert transcript_data.var.loc["tx-known", "gene_id"] == "gene-known"
    assert gene_data.uns["dogme"]["count_type"] == "unique_umi_counts"
    assert gene_data.uns["dogme"]["unmapped_feature_count"] == 1
    assert transcript_data.uns["dogme"]["unmapped_feature_count"] == 1
    assert qc_output.exists()
