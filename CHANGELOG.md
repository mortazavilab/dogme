# Changelog

All notable DOGME changes are documented here.

## [Unreleased]

### Added

- Added optional Dorado barcode demultiplexing through the `kitName` parameter. When configured, DOGME publishes classified BAMs with `${sample}.bcNN` names and keeps unclassified/no-barcode BAMs as publish-only outputs.
- Added barcode-aware mapping, bulk kallisto FASTQ/quantification, RNA/cDNA annotation, and RNA/DNA modkit processing for demultiplexed reads.
- Added `kitName = null` to the bundled configuration examples and test profiles.
- Added optional `keepBarcodes` filtering to retain the specified number of highest-read-count classified barcodes before downstream analysis.
- Added `${sample}.demux_summary.tsv`, reporting read counts, percentages, and downstream-retention status for every Dorado demultiplexing category.

### Fixed

- Ensure the `remap` entry point regenerates SeqSpec and reruns Splitcode for single-cell cDNA inputs from unmapped BAMs.
- Make remap FASTQ extraction consume the discovered BAM input, while preserving support for the existing `${sample}.unmapped.bam` naming convention.
- Correct Parse Evercode WT Mega v2 Nanopore single-cell extraction geometry: recover the 10-base UMI and ordered 24-base barcode from the TruSeq-R2 reverse-complement side, trim poly(T)/TSO technical sequence, and emit normalized cDNA for forward and reverse-complement reads.
- Make splitcode FASTQ combining skip absent orientation streams, exclude ambiguous read IDs before writing, preserve barcode quality orientation, and publish `${sample}_splitcode_qc.tsv` with input, empty, ambiguous, missing-orientation, and emitted-triplet counts.
- Preserve splitcode's configured one-mismatch linker tolerance (`1:1:1`) instead of applying an exact-anchor filter to unclassified reads.
- Demultiplex inline-classified Dorado basecalls with `--no-classify`, preserving classifications after barcode trimming instead of attempting a second classification pass.
- Run QC and inventory reporting only after mapping, quantification, annotation, and modification branches have completed.
- Preserve the configured output directory as the inventory-report input while using completion channels solely for report ordering.
- Suppress undefined-parameter warnings when optional single-cell and SeqSpec configuration values are omitted.
- Render optional SeqSpec command-line arguments before task execution, avoiding invalid literal Groovy expressions in task scripts.
- Add the required `lib_struct` field to bulk cDNA and gDNA SeqSpec templates for SeqSpec 0.3 validation.
- Preserve barcode IDs in annotated BAM, GTF, abundance, QC, and log filenames.
- Preserved barcode identifiers through demultiplexed BAM, FASTQ, SeqSpec, Kallisto, mapping, annotation, and modification outputs.
- Generate a per-barcode FASTQ and SeqSpec artifact for demultiplexed RNA and cDNA reads.
- Allow barcode modification branches with empty mapped strands or no detected modification calls to complete without failing.

### Documentation

- Documented the Parse WT Mega v2 100,000-read empirical regression expectations: 17,308 forward reads, 21,022 reverse-complement reads, 115 ambiguous reads excluded, corrected UMI homopolymer rates near zero, and corrected cDNA-region medians of 463-464 bases.
- Documented the opt-in demultiplexing configuration, output locations, barcode naming, and current single-cell/reporting limitations.

## [1.4.0] - 2026-08-13

### Added

- Added single-cell kallisto/bustools quantification equivalent to the requested `kb count` settings, including long-read mode, custom technology parsing, barcode whitelist correction, and support for precomputed or auto-built `k=63` indexes.

## [1.3.3] - 2026-08-10

### Added

- Optional seqspec generation for FASTQs produced from unmapped BAMs.
- Automatic built-in template selection for RNA, cDNA, DNA, and Parse single-cell configurations.
- Derived FASTQ filename, size, MD5, observed length bounds, and run date, with optional variable overrides.
- Added an explicit Parse Evercode WT mega v2 nanopore seqspec 0.4.0 template with barcode, linker, UMI, and ONT adapter geometry.
- Built-in template rendering followed by `seqspec upgrade`, `seqspec format`, and `seqspec check`.
- Published the final seqspec artifact beside its FASTQ under `${fastqDir}`.
- `singleCell`, `singleCellKit`, `seqspecTemplate`, and `seqspecVariables` configuration parameters.
- Clear runtime errors when the configured Docker/Apptainer image does not provide seqspec.
- Bundled Parse WT v2 and Parse WT Mega v2 barcode-onlist metadata for default single-cell seqspec rendering.
- Single-cell cDNA splitcode processing after seqspec generation, including normalized `ONT.config` geometry, bundled correction/mergeRT passes, and published combined cDNA/UMI/barcode FASTQs under `${fastqDir}/single-cell`.

### Documentation

- Documented seqspec artifact generation from supplied templates and variables. The artifact records declared geometry; it is not inferred from FASTQ or BAM data.
- Documented the Evercode template as an explicit input; its 0.4.0 schema remains unvalidated until it is rendered and checked inside the DOGME image.
- Documented single-cell cDNA splitcode processing, its shared DOGME image runtime requirement, and the splitcode-derived FASTQ handoff for `kb count`.

### Fixed

- FASTQ metadata measurement now supports gzip-compressed `.fastq.gz` inputs.
- Seqspec rendering now supports `seqspec 0.3.1`, including images where `seqspec --version` exits nonzero after reporting its version. Schema 0.3 templates are formatted and checked without the 0.4-only upgrade step.
- Single-cell splitcode correction now supplies paired output paths for cDNA/UMI and barcode processing, avoiding splitcode 0.31.6 output-count failures.
- Single-cell barcode FASTQs now derive from the cDNA correction pass's extracted, on-list barcode sequence before `mergeRT`, keeping the published triplet record-synchronous.
- Single-cell barcode correction uses a generic extraction name instead of the `13G` example label and removes its sample-scoped intermediate at task completion.

## [1.3.2]

### Added

- Standalone `kallisto` entry point for extracting FASTQs from unmapped BAMs.
- cDNA FASTQ input support through the `fastqCDNA` entry point.

### Changed

- Automatic GTF-to-junction BED conversion for RNA and cDNA mapping.
- Automatic kallisto reference construction when prebuilt index and t2g files are not supplied.
- Increased maximum intron size for minimap2 spliced alignment.

### Fixed

- Improved workflow channel grouping and multi-genome handling.

## [1.3.1]

### Changed

- Continued the 1.3.x workflow updates, including two-pass annotation processing and expanded BAM/FASTQ reporting.

### Fixed

- Reduced annotation memory use by processing one chromosome at a time.

## [1.3.0]

### Added

- Standalone reports, remap, basecall, kallisto, and annotation entry points.
- Open-chromatin support for DNA workflows.

### Changed

- Refactored workflow modules and expanded multi-genome processing.
- Added software version reporting and broader output QC.

## [1.2.2]

### Added

- DOGME abundance and GTF outputs from annotation.
- Reconciliation and expanded annotation/QC reporting.

### Changed

- Updated annotation and reconciliation outputs from TALON-style processing to DOGME outputs.

### Fixed

- Corrected novel gene/transcript counts, strand reconciliation, and BAM metadata handling.

## [1.2.1]

### Added

- Basecalling and modification-model selection support.

### Changed

- Automatic Dorado model download and expanded RNA/DNA workflow configuration.
