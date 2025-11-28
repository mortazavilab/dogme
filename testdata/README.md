# Toy Test Data for Dogme
The following folders contain a config file and tiny pod5 files for testing 
Dogme. All of the data is from mouse samples and can be mapped to mm39 or 
equivalent.

Each folder contains a pod5 folder with a few thousand reads in pod5 format.
Each folder also contains two config files. One was tested on the UCI SLURM 
cluster with singularity, the other on a local linux server with docker.

If you will not use either docker or singularity, then install all the 
packages as described by the main dogme readme, and then add all of these
paths in dogme.profile. If using docker/singularity, dogme.profile should be
empty.

In order to run these, you will need to modify the config files to specify:
- Where the genomic fasta and GTF annotations are 
- Where the kallisto-LR index and t2g files are

The genomic fasta and gtf are copies of mm39 from the IGVF Consortium:
https://api.data.igvf.org/reference-files/IGVFFI9282QLXO/@@download/IGVFFI9282QLXO.fasta.gz
https://api.data.igvf.org/reference-files/IGVFFI4777RDZK/@@download/IGVFFI4777RDZK.gtf.gz
