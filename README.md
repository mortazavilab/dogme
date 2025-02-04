# dogme
A nextflow pipeline for basecalling nanopore reads with and without modifications before mapping and processing. Dogme supports 3 different read types: direct RNA ("RNA"), cDNA ("CDNA"), and genomic DNA ("DNA"). Given a folder of pod5 files, Dogme will download the latest dorado models, run dorado, merge bams, and map them with minimap2. For RNA and DNA, it will extract modifications using modkit and filter those. For RNA, cDNA, it will also extract a fastq and run kallisto.
### Installation:
The current version of Dogme has been tested on SLURM clusters and on Macs with the following software versions: 

| Software | Version |
|----------|---------|
| dorado   | 0.9.1   |
| samtools | 1.15.1  |
| minimap2 | 2.28    |
| mod_kit  | 0.4.3   |
| kallisto | 0.51.1  |
| bustools | 0.43.2  |

You can copy and paste this table directly into your GitHub README.md file. The table will be rendered properly in GitHub's markdown viewer.
You must first install dorado, minimap2, samtools, modkit, kallisto, and bustools, and the latest version of nextflow and add their paths to the dogme.profile file in the launch directory or create an empty file of the same name. 
### Prerequisites:
Prerequisites include installing nextflow, using java/17 or better, installing modkit, installing minimap2, installing kallisto (compiled for long-reads), and installing bustools.
By default, nextflow will use the launchdirectory as the place to create its workfolder.
You must also create your own long-read kallisto index and t2g file for your genome of interest. You will need to add the path to the genome fasta, the transcriptome gtf, the kallisto index and the kallisto t2g file to your custom config file. 
### Create your own config file using the template:
The config file must be updated to the list where the top directory containing the pod5 folder and the pipeline outputs are located as well as the work directory. Ideally the config file would also be in the top directory folder. For most use cases only the top parameters in the config file will need to be changed. 
``` 
params {
    sample = 'nxtest'
    //readType can either be 'RNA', 'DNA' or 'CDNA'
    readType = 'RNA'
    // modification type should be set as necessary if different from 'inosine_m6A,pseU,m5C' for RNA and '5mCG_5hmCG,6mA' for DNA. 
    //modifications = 'inosine_m6A,pseU,m5C'
    //change setting if necessary 
    minCov = 3
    perMod = 5
    // change if the launch directory is not where the pod5 and output directories should go
    topDir = "${launchDir}"
    ```
Be sure to change the process section of the example config file to reflect your cluster environment. 
### Running dogme:
Running Dogme on typical dataset can take more than 24 hours, therefore it is recommended to run Dogme within a job or a saved virtual terminal such as screen or Tmux.  Change your folder to be where you want Dogme to run (the 'launch' directory), and launch Dogme directly from github using the following command: 
 ```
  nextflow run mortazavilab/dogme -c yourconfig.conf
  ```
By default, the pipeline will create several folders within the launch directory such as bams, bedMethyl, fastqs, and kallisto - all of which can be customized in the config file. If you need to resume your work add '-resume' to the nextflow command after deleting the html report and trace files.
