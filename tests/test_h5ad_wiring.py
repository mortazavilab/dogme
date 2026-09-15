from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_single_cell_kallisto_writes_gene_and_transcript_h5ad():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    process = nextflow.split("process singleCellKallistoTask", 1)[1].split("process makeKallistoRefsTask", 1)[0]

    assert "params.singleCellH5ad ?" in process
    assert 'bustools count "\\${output_dir}/sorted.bus"' in process
    assert '-o "\\${output_dir}/count" --cm -m -g ${t2gFile}' in process
    assert 'transcript_identity.t2g' in process
    assert '-o "\\${output_dir}/transcript_count" --cm -m \\' in process
    assert '-g "\\${output_dir}/transcript_identity.t2g"' in process
    assert "--feature-type gene" in process
    assert "--feature-type transcript" in process
    assert '${params.sample}_${genomeName}.gene.h5ad' in process
    assert '${params.sample}_${genomeName}.transcript.h5ad' in process
    assert '--t2g "${t2gFile}"' in process


def test_h5ad_parameters_are_defaulted_and_validated():
    nextflow = (ROOT / "dogme.nf").read_text()

    assert "params.singleCellH5ad = params.containsKey('singleCellH5ad') ? params.singleCellH5ad : true" in nextflow
    assert "params.singleCellEntity" in nextflow
    assert "singleCellH5ad must be true or false" in nextflow
    assert "singleCellEntity must be either 'cell' or 'nucleus'" in nextflow


def test_h5ad_runtime_requirements_are_declared():
    requirements = (ROOT / "requirements-h5ad.txt").read_text()

    for package in ("anndata", "h5py", "numpy", "pandas", "scipy"):
        assert package in requirements
