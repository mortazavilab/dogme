from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def test_single_cell_kallisto_writes_gene_and_transcript_h5ad():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    process = nextflow.split("process singleCellKallistoTask", 1)[1].split("process makeKallistoRefsTask", 1)[0]

    assert 'if [[ "${params.singleCellH5ad}" == "true" ]]; then' in process
    assert "path h5adScript" in process
    assert "python ${h5adScript}" in process
    assert 'path "${params.sample}.${genomeName}.kallisto.bam", optional: true, emit: bams' in process
    assert "kallisto bus --long --threshold 0.8" in process
    assert "--genomebam --gtf '${genomeGtf}' --chromosomes chromosomes.txt" in process
    assert 'fasta_index="${genomeFasta}.fai"' in process
    assert 'samtools sort -@ ${task.cpus} -o "${params.sample}.${genomeName}.kallisto.bam"' in process
    assert 'publishDir params.bamDir, pattern: "*.kallisto.bam*", mode: \'copy\', overwrite: true' in process
    assert 'bustools count "\\${output_dir}/sorted.bus"' in process
    assert '-o "\\${output_dir}/count" --cm -m --genecounts -g ${t2gFile}' in process
    assert 'transcript_identity.t2g' in process
    assert '-o "\\${output_dir}/transcript_count" --cm -m --genecounts \\' in process
    assert '-g "\\${output_dir}/transcript_identity.t2g"' in process
    assert "--feature-type gene" in process
    assert "--feature-type transcript" in process
    assert '${params.sample}_${genomeName}.gene.h5ad' in process
    assert '${params.sample}_${genomeName}.transcript.h5ad' in process
    assert '${params.sample}_${genomeName}.single_cell_qc.tsv' in process
    assert "--gene-h5ad" in process
    assert "--transcript-h5ad" in process
    assert "--qc-output" in process
    assert '--t2g "${t2gFile}"' in process
    assert 'file("${projectDir}/scripts/make_h5ad.py")' in nextflow
    assert '.combine(seqspecFile.collect())' in nextflow


def test_h5ad_parameters_are_defaulted_and_validated():
    nextflow = (ROOT / "dogme.nf").read_text()

    assert "def optionalParam(String name, defaultValue = null)" in nextflow
    assert "params.singleCellH5ad = optionalParam('singleCellH5ad', true)" in nextflow
    assert "params.singleCellEntity" in nextflow
    assert "singleCellH5ad must be true or false" in nextflow
    assert "singleCellEntity must be either 'cell' or 'nucleus'" in nextflow


def test_h5ad_runtime_requirements_are_declared():
    requirements = (ROOT / "requirements-h5ad.txt").read_text()

    for package in ("anndata", "h5py", "numpy", "pandas", "scipy"):
        assert package in requirements


def test_seqspec_optional_parameters_are_guarded():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()

    assert "optionalParam('seqspecMd5', true)" in nextflow
    assert "optionalParam('seqspecTemplate')" in nextflow
    assert "optionalParam('seqspecVariables')" in nextflow


def test_genome_bam_generation_is_wired_and_published():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    dogme = (ROOT / "dogme.nf").read_text()

    process = nextflow.split("process genomeBamTask", 1)[1].split("process demuxKallistoTask", 1)[0]
    assert "kallisto quant --long -P ONT --genomebam" in process
    assert "--gtf" in process
    assert "--chromosomes chromosomes.txt" in process
    assert 'path "${sampleName}.${genomeName}.kallisto.bam"' in process
    assert 'path "${sampleName}.${genomeName}.kallisto.bam.bai"' in process
    assert "publishDir params.bamDir" in process
    assert "samtools faidx" in process
    assert 'fasta_index="${genomeFasta}.fai"' in process
    assert 'fasta_index="\\${genomeFasta}.fai"' not in process
    assert "workflow genomeBamWorkflow" in nextflow
    assert "genomeBamReferencesAvailable" in nextflow
    assert "ref.index ?: (refIndex == 0 ? optionalParam('kallistoIndex') : null)" in nextflow
    assert "ref.t2g ?: (refIndex == 0 ? optionalParam('t2g') : null)" in nextflow
    assert "params.genomeBam" in dogme
    assert "!bam.name.endsWith('.kallisto.bam')" in dogme


def test_remap_waits_for_genome_bam_generation():
    nextflow = (ROOT / "nanoporeModule.nf").read_text()
    remap_workflow = nextflow.split("workflow remapWorkflow", 1)[1].split("workflow fastqCDNAWorkflow", 1)[0]
    kallisto_workflow = nextflow.split("workflow kallistoWorkflow", 1)[1].split("workflow kallistoFastqWorkflow", 1)[0]
    genome_bam_workflow = nextflow.split("workflow genomeBamWorkflow", 1)[1].split("// The splitmodification", 1)[0]

    assert "kallistoResults = kallistoWorkflow(unmappedbam)" in remap_workflow
    assert "terminalKallisto.bams.collect()" in kallisto_workflow
    assert "genomeBamWorkflow(fastqFile.map { fastq -> tuple(fastq, params.sample) }).completion" in kallisto_workflow
    assert "def terminalKallistoCompletion" in kallisto_workflow
    assert "genomeBamCompletion = genomeBamTask.out.bams.collect()" in genome_bam_workflow
    assert "completion = genomeBamCompletion" in genome_bam_workflow
