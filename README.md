# dogme

A nextflow pipeline for basecalling nanopore reads with and without modifications before mapping and processing. Dogme supports 3 different read types: direct RNA ("RNA"), cDNA ("CDNA"), and genomic DNA ("DNA"). Given a folder of pod5 files, Dogme will download the latest dorado models, run dorado, merge bams, and map them with minimap2. For RNA and DNA, it will extract modifications using modkit and filter those. For RNA and cDNA, it will also extract a fastq and run kallisto.

---

## What's New in Dogme 1.4.0

- **End-to-end single-cell cDNA workflow:** With `readType = 'CDNA'` and `singleCell = true`, DOGME generates and validates seqspec metadata, splits and corrects cDNA/UMI/barcode FASTQs with splitcode, and runs single-cell kallisto/bustools quantification.
- **FASTQ seqspec generation:** DOGME can render, upgrade, format, and validate a seqspec artifact whenever it generates a FASTQ from an unmapped BAM. Single-cell cDNA runs additionally generate a splitcode configuration and processed FASTQ.
- **Optional Dorado demultiplexing:** Set `kitName` to classify and demultiplex barcoded reads. Classified reads use `${sample}.${barcode}` output names and run through barcode-aware mapping, bulk kallisto, annotation, and RNA/DNA modification processing; unclassified reads are published without downstream analysis.

- **fastq CDNA support**: 
  New workflow and entry point to start from an existing CDNA FASTQ, create an unmapped BAM, run minimap2 and kallisto in parallel, then annotate the mapped BAMs.
- **Modkit Open Chromatin support:**  
  New workflow and entry point to call open chromatin signal and regions in mapped BAMs with m6A modifications using `modkit 0.5` and higher. This produces both a bed files of regions and a bedgraph per genome.
- **Transcript Annotation:**  
  New workflow entry point to annotate mapped BAMs with transcript information using `annotateRNA.py`. This produces annotated BAM files, TALON outputs and QC summary CSVs for each genome.
- **Standalone kallisto entry point:**  
  Added a `kallisto` entry point that starts from unmapped BAMs, extracts FASTQ, and runs the kallisto long-read quantification steps without re-running basecalling or remapping.
- **Single-cell kallisto/bustools quantification:**
  Added a `kb count`-equivalent workflow using kallisto and bustools alone, with long-read mode, technology string `2,0,24:1,0,10:0,0,0`, barcode whitelist correction, and support for precomputed or auto-built `k=63` indexes.
- **Automatic GTF-to-Junction BED Conversion:**  
  The pipeline now automatically converts GTF files to junction BED files for minimap2 spliced alignment, ensuring correct handling of RNA and cDNA mapping.
- **Increased Maximum Intron Size:**  
  The minimap2 mapping step now uses `--splice-max 500000` for improved detection of long introns in spliced alignments.
- **Improved Workflow Modularity:**  
  New entry points / workflows to improve modularity and restartability: basecall, remap, reports, kallisto (start-from-unmapped-BAMs), annotateRNA (start-from-mapped-BAMs), and the original main workflow.
- **Reporting and QC:**
  - generate_report.py now gathers additional per-BAM and per-FASTQ statistics into qc_summary.csv and inventory_report.tsv.
  - Reports workflow allows generation of metadata/QC without re-running basecalling or mapping.
- **Bug Fixes and Robustness:**  
  - Lower memory footprint for annotateRNA.py using 2-pass processing first to count then to write the bam file on a per-chromosome basis.
  - Improved handling of file naming, channel grouping, and tuple passing so multi-genome and multi-strand runs behave correctly.
  - dorado model download is only run if the model directory does not already exist (avoids repeated downloads).
  - Processes include retry/error strategies for robustness of long-running tasks.


Dogme 1.4.0 carries forward the 1.3.3 workflow updates, including seqspec generation, single-cell cDNA splitting, and single-cell kallisto/bustools quantification.

---

## Dogme helper python script

The following Python scripts are included or updated in the scripts/ directory. These are referenced by the Nextflow processes and are important to the new features:

- scripts/software_versions.py
  - Collects and records the software and model versions used for a run into ${sample}.softwareVersion.txt.

- scripts/gtf_to_junction_bed.py
  - Converts a GTF into a junction BED suitable for minimap2 spliced alignment. This is used automatically for RNA and CDNA mapping.

- scripts/filterbed.py
  - Filters modkit bed outputs by minimum coverage and per-mod thresholds (params.minCov and params.perMod), reading and writing compressed `.bed.gz` files for the BED workflow.

- scripts/annotateRNA.py
  - Annotates mapped BAMs with transcript information. Now outputs TALON-compatible outputs and expanded QC CSVs by default. Accepts a -CDNA flag for cDNA-specific behavior.

- scripts/generate_report.py
  - Gathers inventory and QC metrics across outputs and generates inventory_report.tsv and qc_summary.csv (additional BAM/FASTQ stats added post-1.2.2), including compressed `bedMethyl` outputs.

- scripts/reconcileBams.py
  - Consolidates per-sample/TALON-style BAM/gene outputs and reports consolidated gene counts and reconciliation statistics (updated in 1.2.2 to correct gene counts).

- scripts/report_bed_rep.py
  - Compares per-sample modkit BED outputs across folders (replicates), including compressed `.bed.gz` files, merges plus/minus strands per sample and produces per-reference/modification comparison CSVs (e.g., Compare_GRCh38_m6A.csv).

---

## Installation

The current version of Dogme has been tested on SLURM clusters and on Macs with the following software versions: 

| Software | Version |
|----------|---------|
| dorado   | 1.2     |
| samtools | 1.15.1  |
| minimap2 | 2.28    |
| mod_kit  | 0.5     |
| kallisto | 0.51.1  |
| bustools | 0.43.2  |

You must first install dorado, minimap2, samtools, modkit, kallisto, and bustools, and the latest version of nextflow and add their paths to the dogme.profile file in the launch directory or create an empty file of the same name. 

Alternatively, you can use the provided Docker/Singularity image as described in the Container image section below.

---

## Prerequisites

Prerequisites include installing nextflow, using java/17 or better, installing modkit, installing minimap2, installing kallisto (compiled for long-reads), and installing bustools.
By default, nextflow will use the launchdirectory as the place to create its workfolder.
You can either provide pre-built kallisto index and t2g files in your config, or set `kallistoIndex` and `t2g` to `null` to have the pipeline auto-build them from the genome and GTF files in your `genome_annot_refs`. When auto-building, the pipeline will generate a kallisto index and t2g for each genome entry and run quantification against each one.

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

    // Optional seqspec generation for FASTQs created by DOGME.
    singleCell = false
    singleCellKit = null
    seqspecTemplate = null
    seqspecVariables = null
    // Optional Dorado barcode kit identifier. Leave null to disable demultiplexing.
    kitName = null

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
    // Or set both to null to auto-build from genome_annot_refs:
    // kallistoIndex = null
    // t2g = null
    
    //default accuracy is sup
    accuracy = "sup"
    // change this value if 0.9 is too strict
    // if set to null or '' then modkit will determine its threshold by sampling reads. 
    modkitFilterThreshold = 0.9

    //rest of config file - see dogmetest-param.conf
```
  
Be sure to change the process section of the example config file to reflect your cluster environment. 

### Optional Dorado demultiplexing

Set `kitName` to the Dorado kit identifier to enable barcode classification and demultiplexing during the `main` or `basecall` workflow:

```
kitName = 'SQK-RNA004'
```

Leave `kitName` unset, `null`, or blank to preserve the normal single-sample basecalling behavior. When enabled, Dorado publishes barcode BAMs below `${dorDir}/demux` and classified BAMs below `${bamDir}` using the `${sample}.${barcode}` prefix. Classified RNA and cDNA reads are mapped, extracted to barcode-specific FASTQs, quantified with bulk kallisto, and annotated. Classified RNA and DNA reads also run through barcode-specific modkit pileup and BED filtering. Unclassified and no-barcode BAMs are published but are not analyzed. Single-cell cDNA routing and run-level aggregate reports are not yet enabled for demultiplexed runs. The exact kit identifier must be supported by the Dorado version installed in the execution container or profile.

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
By default, the pipeline will create several folders within the launch directory such as bams, bedMethyl, fastqs, and kallisto - all of which can be customized in the config file. The `bedMethyl` directory stores compressed `.bed.gz` outputs for the raw modkit BEDs, filtered BEDs, and final per-modification BEDs. If you need to resume your work add '-resume' to the nextflow command after deleting the html report and trace files.

---

## Entry Points

- **main**: Full pipeline from pod5 files to mapped/annotated/filtered outputs.
- **basecall**: basecall pod5 files into unmapped bam file.
- **remap**: Remap reads starting from unmapped BAM files.
- **fastqCDNA**: Start from an existing CDNA FASTQ, create an unmapped BAM, run minimap2 and kallisto in parallel, then annotate the mapped BAMs.
- **kallisto**: Extract FASTQ from unmapped BAM files and run the kallisto long-read quantification steps.
- **modkit**: Run modification extraction and filtering.
- **reports**: Generate summary reports only.
- **annotateRNA**: Annotate mapped BAMs with transcript information and produce QC summaries.

## Seqspec generation for generated FASTQs

When DOGME creates a FASTQ from an unmapped BAM, it automatically creates a seqspec artifact beside the FASTQ. Region structure is declared by a built-in assay template; read lengths and file metadata are measured from the generated FASTQ. An unmapped BAM carries no assay geometry, so the artifact does not infer region structure from the BAM.

Single-cell read processing is enabled for `readType = 'CDNA'` with `singleCell = true`. After seqspec generation, DOGME creates `ONT.config`, normalizes splitcode's `3:3:3` geometry to `1:1:1`, and runs `splitcode` with two threads. The original FASTQ remains available, and the combined `${sample}_barcode.fastq.gz`, `${sample}_cDNA.fastq.gz`, and `${sample}_umi.fastq.gz` files are published under `${fastqDir}/single-cell`.

The task combines the splitcode `c_*`, `f_*`, `r_*`, and `rc_*` cDNA, UMI, and barcode files, then applies the bundled `config-correct.txt` and `config.mergeRT` splitcode passes. The cDNA correction pass extracts the recognized on-list barcode segments, and that sample-scoped extracted FASTQ is passed to `mergeRT` before publication. It publishes three final FASTQs: `${sample}_barcode.fastq.gz`, `${sample}_cDNA.fastq.gz`, and `${sample}_umi.fastq.gz`. With `singleCell = true`, the kallisto workflow consumes these files directly and expands the equivalent of `kb count` into kallisto and bustools commands using technology string `2,0,24:1,0,10:0,0,0`, `--long`, and `--threshold 0.8`. Barcode correction and `CB`/`CR`/`UB`/`UR` tagging are not performed by this splitcode task; bustools whitelist correction is performed by the single-cell quantification task.

The built-in template is selected from the run configuration:

| `readType` | `singleCell` | Template |
| --- | --- | --- |
| `RNA` | `false` | `templates/seqspec/ont-bulk-drna.yaml.j2` |
| `CDNA` | `false` | `templates/seqspec/ont-bulk-cdna.yaml.j2` |
| `DNA` | `false` | `templates/seqspec/ont-bulk-gdna.yaml.j2` |
| `CDNA` | `true` | `templates/parse-evercode-wt-mega-v2-nanopore.yaml.j2` |

The normal case requires no seqspec-specific parameters. An unsupported
`readType`/`singleCell` combination fails clearly rather than using the wrong
assay geometry. `params.seqspecTemplate` remains an explicit override, and
`params.seqspecVariables` is an optional JSON object merged over values derived
from the run.

The derived values include the generated FASTQ filename, file ID, byte size,
local URL, MD5 checksum (unless `params.seqspecMd5 = false`), observed minimum
and maximum read lengths, and run date. User-provided variables take precedence
over those derived values.

Explicit override configuration remains available:

```groovy
params {
  singleCell = false
  singleCellKit = null
  seqspecTemplate = '/path/to/assay-seqspec.yaml.j2'
  seqspecVariables = '/path/to/assay-seqspec-variables.json'
}
```

For built-in Parse single-cell specs, `singleCellKit` defaults to `parse-wt-mega-v2`. Set it to either `parse-wt-v2` or `parse-wt-mega-v2` to populate the bundled barcode onlist URLs, checksums, assay metadata, and ONT sequencing metadata. User-provided `seqspecVariables` take precedence over those defaults. The remote onlist file sizes are unavailable in the bundled metadata and are rendered as `null` unless supplied in `seqspecVariables`.

The repository includes the Parse Evercode WT mega v2 nanopore template at
`templates/parse-evercode-wt-mega-v2-nanopore.yaml.j2`. It is selected for
`readType = 'CDNA'` with `singleCell = true`. Its 0.4.0 schema has not yet been
validated inside the DOGME image; run `seqspec check` after rendering before
using it for production data.

The renderer and splitcode task run inside `ghcr.io/mortazavilab/dogme-pipeline:latest` by default. The image must provide both `seqspec` and `splitcode` on `PATH`. DOGME fills template placeholders without Jinja2. Enable Docker or Singularity/Apptainer in the Nextflow configuration; without a container runtime, the tasks cannot access these image-provided dependencies. Seqspec artifacts are published beside the generated FASTQ under `${fastqDir}`, while splitcode outputs are published under `${fastqDir}/single-cell`.

`singleCell` defaults to `false`, and splitcode runs only for `readType = 'CDNA'` with `singleCell = true`. A pre-rendered external spec may be supplied with `params.seqspec` for workflows that consume existing FASTQs.

---

## Example: Annotating BAMs

To annotate mapped BAMs with transcript information:

```
nextflow run mortazavilab/dogme -entry annotateRNA -c yourconfig.conf
```

This will produce annotated BAMs and QC summary files for each genome using the mapped bams.

## Example: Running kallisto From Unmapped BAMs

To run the extracted FASTQ plus kallisto quantification steps starting from existing unmapped BAMs:

```
nextflow run mortazavilab/dogme -entry kallisto -c yourconfig.conf
```

This will look for `*.unmapped.bam` files in `params.bamDir` and write kallisto outputs without re-running basecalling or remapping.

## Example: Running fastqCDNA From an Existing CDNA FASTQ

To start from a pre-existing CDNA FASTQ, place exactly one file named `${params.sample}.fastq.gz` or `${params.sample}.fastq` in `params.fastqDir` and run:

```
nextflow run mortazavilab/dogme -entry fastqCDNA -c yourconfig.conf
```

This entry creates `${params.sample}.unmapped.bam` in `params.bamDir`, runs kallisto in parallel with minimap2, and then annotates the mapped BAM outputs.

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
