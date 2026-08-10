# Changelog

All notable DOGME changes are documented here.

## [1.3.3] - 2026-08-10

### Added

- Optional seqspec generation for FASTQs produced from unmapped BAMs.
- Jinja2 template rendering followed by `seqspec upgrade`, `seqspec format`, and `seqspec check`.
- Published seqspec artifacts and render variables under `${fastqDir}/seqspec`.
- `singleCell`, `singleCellKit`, `seqspecTemplate`, and `seqspecVariables` configuration parameters.
- Explicit single-cell opt-in: only the boolean `singleCell = true` enables single-cell behavior.
- Clear runtime errors when the configured Docker/Apptainer image does not provide seqspec or Jinja2.

### Documentation

- Documented seqspec generation for bulk and single-cell FASTQ workflows.
- Documented the separation between seqspec assay geometry and assay-specific barcode correction.

## [1.3.2]

- Previous DOGME release.