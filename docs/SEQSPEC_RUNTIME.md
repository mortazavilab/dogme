# Seqspec Runtime

DOGME invokes seqspec inside the configured Docker/Apptainer image. Template filling is implemented by DOGME and does not require Jinja2. The image recipe is maintained outside this repository; its authoritative location is **TODO: record the external image recipe URL or repository**.

The seqspec version shipped by the current image is **TODO: fill from `seqspec --version` during item 1a validation**.

The renderer must run in an image that provides `seqspec`, `seqspec upgrade`, `seqspec format`, and `seqspec check`. Host installation of seqspec is not required for containerized DOGME runs.