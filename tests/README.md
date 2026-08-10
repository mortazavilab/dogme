# Tests

Run the test harness with:

```bash
pytest -q
```

`test_generate_seqspec.py` skips when `seqspec` or Jinja2 is unavailable. Run it inside the DOGME Docker/Apptainer image after item 1a validates the image runtime.

The fixture `tests/fixtures/synthetic-seqspec.yaml.j2` has never been through `seqspec check`. If item 1a fails, it will therefore be ambiguous whether the failure comes from the renderer or from the fixture itself.