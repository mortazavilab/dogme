# Changelog

All notable DOGME changes are documented here.

## [1.3.3] - 2026-08-10

### Added

- Optional seqspec generation for FASTQs produced from unmapped BAMs.
- Added an explicit Parse Evercode WT mega v2 nanopore seqspec 0.4.0 template with barcode, linker, UMI, and ONT adapter geometry.
- Jinja2 template rendering followed by `seqspec upgrade`, `seqspec format`, and `seqspec check`.
- Published seqspec artifacts and render variables under `${fastqDir}/seqspec`.
- `singleCell`, `singleCellKit`, `seqspecTemplate`, and `seqspecVariables` configuration parameters.
- Clear runtime errors when the configured Docker/Apptainer image does not provide seqspec or Jinja2.

### Documentation

- Documented seqspec artifact generation from supplied templates and variables. The artifact records declared geometry; it is not inferred from FASTQ or BAM data.
- Documented the Evercode template as an explicit input; its 0.4.0 schema remains unvalidated until it is rendered and checked inside the DOGME image.
- Documented that single-cell read processing, barcode/UMI extraction, `CB`/`CR`/`UB`/`UR` tagging, and barcode correction are not implemented. `singleCell` and `singleCellKit` remain reserved parameters.

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
