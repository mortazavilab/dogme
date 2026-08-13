# Changelog

All notable DOGME changes are documented here.

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
- Single-cell splitcode correction now supplies the required two output paths for paired FASTQ processing while selecting only the primary output, avoiding splitcode 0.31.6 output-count failures.

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
