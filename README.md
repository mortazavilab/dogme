# dogme

A nextflow pipeline for basecalling nanopore reads with and without modifications before mapping and processing. Dogme supports 3 different read types: direct RNA ("RNA"), cDNA ("CDNA"), and genomic DNA ("DNA"). Given a folder of pod5 files, Dogme will download the latest dorado models, run dorado, merge bams, and map them with minimap2. For RNA and DNA, it will extract modifications using modkit and filter those. For RNA and cDNA, it will also extract a fastq and run kallisto.

---

## What's New in Dogme 1.2

- **Modkit Open Chromatin support:**  
  New workflow and entry point to call open chromatin signal and regions in mapped BAMs with m6A modifications using `modkit 0.5` and higher. This produces both a bed files of regions and a bedgraph per genome.
- **Transcript Annotation:**  
  New workflow entry point to annotate mapped BAMs with transcript information using `annotateRNA.py`. This produces annotated BAM files, TALON outputs and QC summary CSVs for each genome.
- **Automatic GTF-to-Junction BED Conversion:**  
  The pipeline now automatically converts GTF files to junction BED files for minimap2 spliced alignment, ensuring correct handling of RNA and cDNA mapping.
- **Increased Maximum Intron Size:**  
  The minimap2 mapping step now uses `--splice-max 500000` for improved detection of long introns in spliced alignments.
- **Improved Workflow Modularity:**  
  New entry points and modular workflows allow starting from mapped BAMs, remapping, or running only annotation/reporting as needed.
- **Bug Fixes and Robustness:**  
  Improved handling of file naming, tuple passing, and channel joining for multi-genome and multi-strand workflows.

Dogme 1.2.2 updates `annotateRNA.py` to output the old optional TALON outputs as dogme defaults, updated `reconcileBams.py` to report the correct number of consolidated genes, and added additional statistics for bams and fastqs in the final qc summary. 

---

## What's New in 1.2.3

The pipeline has continued to evolve after 1.2. Notable additions and changes in the post-1.2 releases:

- Workflow / mapping
  - Added a dedicated gtf -> junction BED conversion process (gtf_to_junction_bed.py) and automatic use of junction BEDs for spliced minimap2 alignment.
  - Minimap2 spliced mapping uses a larger maximum intron size (`--splice-max 500000`) to better detect very long introns.
  - New entry points / workflows to improve modularity and restartability: basecall, remap, reports, annotateRNA (start-from-mapped-BAMs), and the original main workflow.
  - Improved grouping and handling of multi-genome runs: mapped BAMs are grouped by genome before annotation to ensure correct pairing with genome GTFs.

- Modkit / open chromatin
  - Per-chromosome open-chromatin calling for DNA (modkit 0.5+) with both per-chromosome bed and bedgraph outputs and downstream consolidation steps.
  - Consolidation tasks combine per-chromosome bed/bg outputs into per-genome files:
    - ${sample}.${genome}.m6Aopen.bed
    - ${sample}.${genome}.m6Aopen.bg
  - Modkit pileup call in the pipeline now accepts an optional threshold parameter (params.modkitFilterThreshold) and is used automatically when provided.

- Reporting and QC
  - generate_report.py now gathers additional per-BAM and per-FASTQ statistics into qc_summary.csv and inventory_report.tsv.
  - Reports workflow allows generation of metadata/QC without re-running basecalling or mapping.

- Annotation improvements
  - annotateRNA.py remains the annotation engine; the pipeline now groups BAMs by genome and pairs them with the correct GTF before invoking the annotator.
  - annotateRNATask supports a CDNA option and produces annotated BAMs, TALON outputs and per-genome QC CSVs by default.

- Robustness / misc
  - Improved handling of file naming, channel grouping, and tuple passing so multi-genome and multi-strand runs behave correctly.
  - dorado model download is only run if the model directory does not already exist (avoids repeated downloads).
  - Processes include retry/error strategies for robustness of long-running tasks.

---

## Dogme helper python script

The following Python scripts are included or updated in the scripts/ directory. These are referenced by the Nextflow processes and are important to the new features:

- scripts/software_versions.py
  - Collects and records the software and model versions used for a run into ${sample}.softwareVersion.txt.

- scripts/gtf_to_junction_bed.py
  - Converts a GTF into a junction BED suitable for minimap2 spliced alignment. This is used automatically for RNA and CDNA mapping.

- scripts/filterbed.py
  - Filters modkit bed outputs by minimum coverage and per-mod thresholds (params.minCov and params.perMod).

- scripts/annotateRNA.py
  - Annotates mapped BAMs with transcript information. Now outputs TALON-compatible outputs and expanded QC CSVs by default. Accepts a -CDNA flag for cDNA-specific behavior.

- scripts/generate_report.py
  - Gathers inventory and QC metrics across outputs and generates inventory_report.tsv and qc_summary.csv (additional BAM/FASTQ stats added post-1.2.2).

- scripts/reconcileBams.py
  - Consolidates per-sample/TALON-style BAM/gene outputs and reports consolidated gene counts and reconciliation statistics (updated in 1.2.2 to correct gene counts).

- scripts/report_bed_rep.py
  - Compares per-sample modkit BED outputs across folders (replicates), merges plus/minus strands per sample and produces per-reference/modification comparison CSVs (e.g., Compare_GRCh38_m6A.csv).

---

## Installation

The current version of Dogme has been tested on SLURM clusters and on Macs with the following software versions: 

| Software | Version |
|----------|---------|
| dorado   | 1.0     |
| samtools | 1.15.1  |
| minimap2 | 2.28    |
| mod_kit  | 0.5     |
| kallisto | 0.51.1  |
| bustools | 0.43.2  |

You must first install dorado, minimap2, samtools, modkit, kallisto, and bustools, and the latest version of nextflow and add their paths to the dogme.profile file in the launch directory or create an empty file of the same name. 

---

## Prerequisites

Prerequisites include installing nextflow, using java/17 or better, installing modkit, installing minimap2, installing kallisto (compiled for long-reads), and installing bustools.
By default, nextflow will use the launchdirectory as the place to create its workfolder.
You must also create your own long-read kallisto index and t2g file for your genome of interest. You will need to add the path to the genome fasta, the transcriptome gtf, the kallisto index and the kallisto t2g file to your custom config file. 

---

## Create your own config file using the template

The config file must be updated to the list where the top directory containing the pod5 folder and the pipeline outputs are located as well as the work directory. Ideally the config file would also be in the top directory folder. For most use cases only the top parameters in the config file will need to be changed. 

``` 
params {
    sample = 'C2C12_RNA_r1'
    //readType can either be 'RNA', 'DNA' or 'CDNA'
    readType = 'RNA'
    // modification type should be set as necessary if different from 'inosine_m6A,pseU,m5C' for RNA and '5mCG_5hmCG,6mA' for DNA. 
    //modifications = 'inosine_m6A,pseU,m5C'
    //change setting if necessary 
    minCov = 3
    perMod = 5
    // change if the launch directory is not where the pod5 and output directories should go
    topDir = "${launchDir}"

    // the following file should be edited to add all the necessary paths for commands such as
    // dorado, samtools, minimap2, kallisto, and bustools
    scriptEnv = "${launchDir}/dogme.profile"

    // needs to be modified to match the right genomic reference
     genome_annot_refs = [
     [name: 'mm39', genome: '/path/genomeRef/IGVFFI9282QLXO.fasta', annot: '/path/genomeRef/IGVFFI4777RDZK.gtf'],
     [name: 'C57BL_6J_T2T_v1', genome: '/path/genomeRef/C57BL_6J_T2T_v1/unmasked.fa', annot: '/path/genomeRef/C57BL_6J_T2T_v1/genes.gtf'],
     [name: 'CAST_EiJ_T2T_v1', genome: '/path/genomeRef/CAST_EiJ_T2T_v1/unmasked.fa', annot: '/path/genomeRef/CAST_EiJ_T2T_v1/genes.gtf']
     ]
    kallistoIndex = '/path/kallistoref/mm39GencM36_k63.idx'
    t2g = '/pathA/kallistoref/mm39GencM36_k63.t2g'
    
    //default accuracy is sup
    accuracy = "sup"
    // change this value if 0.9 is too strict
    // if set to null or '' then modkit will determine its threshold by sampling reads. 
    modkitFilterThreshold = 0.9

    //rest of config file - see dogmetest-param.conf
```
  
Be sure to change the process section of the example config file to reflect your cluster environment. 

---

## Running modkit to call open chromatin in DNA

Dogme will use modkit 0.5+ with the delivered models in GPU mode using the GPU library, based on the version in the path and several shell script variables defined in your local dogme.profile:
```
export MODKITBASE=/path/to/modkit
export MODKITMODEL=${MODKITBASE}/models/r1041_e82_400bps_hac_v5.2.0@v0.1.0
# the following are only needed if using the torch version of modkit
export LIBTORCH=${MODKITBASE}/libtorch
export DYLD_LIBRARY_PATH=${LIBTORCH}/lib
export LD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}
export PATH=${MODKITBASE}/dist_modkit_v0.5.0_5120ef7_tch:$PATH
```

---

## Running dogme

Running Dogme on typical dataset can take more than 24 hours, therefore it is recommended to run Dogme within a job or a saved virtual terminal such as screen or Tmux.  Change your folder to be where you want Dogme to run (the 'launch' directory), and launch Dogme directly from github using the following command: 

 ```
nextflow run mortazavilab/dogme -c yourconfig.conf
```
By default, the pipeline will create several folders within the launch directory such as bams, bedMethyl, fastqs, and kallisto - all of which can be customized in the config file. If you need to resume your work add '-resume' to the nextflow command after deleting the html report and trace files.

---

## Entry Points

- **main**: Full pipeline from pod5 files to mapped/annotated/filtered outputs.
- **basecall**: basecall pod5 files into unmapped bam file.
- **remap**: Remap reads starting from unmapped BAM files.
- **modkit**: Run modification extraction and filtering.
- **reports**: Generate summary reports only.
- **annotateRNA**: Annotate mapped BAMs with transcript information and produce QC summaries.

---

## Example: Annotating BAMs

To annotate mapped BAMs with transcript information:

```
nextflow run mortazavilab/dogme -entry annotateRNA -c yourconfig.conf
```

This will produce annotated BAMs and QC summary files for each genome using the mapped bams.

---

## Container support (Singularity / Docker)

Dogme includes a published container image to simplify reproducible runs and avoid installing all dependencies by hand:

- Container image: ghcr.io/mortazavilab/dogme-pipeline:latest

Two common modes:

- Singularity / Apptainer (recommended on HPC)
  - Enable in your Nextflow config: `singularity.enabled = true`
  - Specify the container and bind mounts in `process` config (example shown in SingularityConfigExample.conf).
  - GPU access: add `--nv` to `process.containerOptions` for GPU-enabled steps (dorado).
  - Example snippet (from SingularityConfigExample.conf):
    ```groovy
    singularity {
        enabled = true
        autoMounts = true
    }

    process {
        container = 'ghcr.io/mortazavilab/dogme-pipeline:latest'
        containerOptions = "--bind /path/to/your/data1,/path/to/your/data2"
        beforeScript = 'export PATH=/opt/conda/bin:$PATH'

        withName: 'doradoTask' {
            containerOptions = "--nv --bind /path/to/your/data1,/path/to/your/data2"
        }
    }
    ```

- Docker (local/workstation)
  - If using Docker, disable Singularity in nextflow and set `process.container` to the same image. Use Docker runtime options via `process.containerOptions`, e.g. `--gpus all` for GPU support.
  - Example:
    ```
    process {
        container = 'ghcr.io/mortazavilab/dogme-pipeline:latest'
        containerOptions = "--gpus all -v /path/to/your/data1:/path/to/your/data1 -v /path/to/your/data2:/path/to/your/data2"
    }
    ```
  - Note: on Macs Docker behaves differently for GPU — GPUs are typically not available on macOS Docker; use a Linux/GPU host or Singularity on cluster for GPU tasks.

Important notes:
- If you run Dogme inside the container image above you do not need to install the listed tools on the host. If you choose not to use containers, you must install dorado, minimap2, samtools, modkit, kallisto, and bustools and ensure they are visible in the PATH (see dogme.profile).
- Ensure any large shared storage mountpoints used by the pipeline (e.g. /path/to/your/data1, /path/to/your/data2) are bound into the container with `containerOptions` so the container can read/write data.
- GPU-enabled steps (dorado / modkit GPU) require adding `--nv` for Singularity or `--gpus` for Docker and a host with GPUs + appropriate drivers.
