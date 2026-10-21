from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_single_cell_kallisto_keeps_bus_quantification_without_genome_bam():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    process = nextflow.split("process singleCellKallistoTask", 1)[1].split("process makeKallistoRefsTask", 1)[0]

    assert "kallisto bus --long --threshold 0.8" in process
    assert '-o "\\${output_dir}/count" --cm -m --genecounts -g ${t2gFile}' in process
    assert "--genomebam" not in process
    assert ".kallisto.bam" not in process
    assert "h5adScript" not in process


def test_annotated_bam_tags_drive_single_cell_h5ad():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    tag_process = nextflow.split("process singleCellAnnotatedBamTask", 1)[1].split("workflow modificationWorkflow", 1)[0]
    annotation_workflow = nextflow.split("workflow annotateRNAWorkflow", 1)[1].split("workflow modificationWorkflow", 1)[0]

    assert "--umi" in tag_process
    assert "--barcode" in tag_process
    assert '--bam "\\${tagged_bam}"' in tag_process
    assert '--gtf "${dogmeGtf}"' in tag_process
    assert "--gene-output" in tag_process
    assert "--transcript-output" in tag_process
    assert 'file("${projectDir}/scripts/tag_annotated_bam.py")' in annotation_workflow
    assert "annotationResults.out.annotation_products" in annotation_workflow
    assert "singleCellAnnotatedBamTask(" in annotation_workflow
    assert "kallistoResults.fastqs" in nextflow
    assert '"${sampleName}.unmapped.${genomeName}.annotated.bam"' in annotation_workflow


def test_h5ad_parameters_are_defaulted_and_validated():
    nextflow = (ROOT / "dogme.nf").read_text()

    assert "def optionalParam(String name, defaultValue = null)" in nextflow
    assert "params.singleCellH5ad = optionalParam('singleCellH5ad', true)" in nextflow
    assert "params.singleCellEntity" in nextflow
    assert "singleCellH5ad must be true or false" in nextflow
    assert "singleCellEntity must be either 'cell' or 'nucleus'" in nextflow
    assert "genomeBam" not in nextflow


def test_h5ad_runtime_requirements_are_declared():
    requirements = (ROOT / "requirements-h5ad.txt").read_text()

    for package in ("anndata", "h5py", "numpy", "pandas", "scipy", "pysam"):
        assert package in requirements


def test_genome_bam_generation_is_removed_but_stale_outputs_are_filtered():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    dogme = (ROOT / "dogme.nf").read_text()

    assert "process genomeBamTask" not in nextflow
    assert "workflow genomeBamWorkflow" not in nextflow
    assert "genomeBamReferencesAvailable" not in nextflow
    assert "--genomebam" not in nextflow
    assert "genomeBam" not in dogme
    assert "!bam.name.endsWith('.kallisto.bam')" in dogme


def test_remap_and_fastq_workflows_connect_splitcode_to_annotation():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    remap_workflow = nextflow.split("workflow remapWorkflow", 1)[1].split("workflow fastqCDNAWorkflow", 1)[0]
    fastq_workflow = nextflow.split("workflow fastqCDNAWorkflow", 1)[1].split("workflow reportsWorkflow", 1)[0]
    kallisto_workflow = nextflow.split("workflow kallistoWorkflow", 1)[1].split("workflow kallistoFastqWorkflow", 1)[0]

    assert "kallistoResults = kallistoWorkflow(unmappedbam)" in remap_workflow
    assert "kallistoResults.fastqs" in remap_workflow
    assert "kallistoResults.fastqs" in fastq_workflow
    assert "singleCellTaggingEnabled" in nextflow
    assert "splitcodeFastqs = splitcodeTask.out.fastqs.map" in kallisto_workflow
    assert "def genomeName = params.genome_annot_refs[0].name" in kallisto_workflow
    assert "completion = terminalKallistoCompletion.combine(seqspecFile.collect())" in kallisto_workflow


def test_seqspec_optional_parameters_are_guarded():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()

    assert "optionalParam('seqspecMd5', true)" in nextflow
    assert "optionalParam('seqspecTemplate')" in nextflow
    assert "optionalParam('seqspecVariables')" in nextflow
