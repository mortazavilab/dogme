# Tests

Run the test harness with:

```bash
pytest -q
```

`test_generate_seqspec.py` skips only the integration test when `seqspec` is unavailable. Template filling is implemented by DOGME and does not require Jinja2. Run the integration test inside the DOGME Docker/Apptainer image after item 1a validates the image runtime.

The fixture `tests/fixtures/synthetic-seqspec.yaml.j2` has never been through `seqspec check`. If item 1a fails, it will therefore be ambiguous whether the failure comes from the renderer or from the fixture itself.

The Parse Evercode 0.4.0 template at `templates/parse-evercode-wt-mega-v2-nanopore.yaml.j2` is also unvalidated until it has been rendered and checked inside the DOGME image. Its barcode on-list variables must refer to real, accessible resources for an image validation run.

The H5AD and BAM-tagging unit tests require the packages in `requirements-h5ad.txt`. Install that file in the test environment or run the tests inside an image containing those packages. The dependency-independent H5AD wiring tests run without them. Single-cell tests verify read-ID matching, `CB`/`UB` tags, UMI deduplication, and inclusion of novel gene/transcript IDs from annotated BAM tags.