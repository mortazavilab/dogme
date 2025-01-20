# dogme
A nextflow pipeline for modification calling in nanopore reads. Given a folder of pod5 files, it will download dorado models, run dorado, merge bams, map them with minimap, extract modifications using modkit and filter those. For RNA, it will also extract a fastq and run kallisto.
### Installation:
The current version of Dogme is customized to run on SLURM clusters. You must first install Dorado, Minimap2, samtools, modkit, kallisto, and bustools and add their paths to the dogme.profile file in the launch directory or create an empty file of the same name. 
### Prerequisites:
Prerequisites include installing nextflow, using java/17 or better, installing dogme in an accessible folder, installing modkit, installing minimap2, installing kallisto (compiled for long-reads), and install bustools.
By default, nextflow will use the launchdirectory as the place to create its workfolder.
You must also create your own long-read kallisto index and t2g file for your genome of interest. You will need to add the path to the genome fasta, the transcriptome gtf, the kallisto index and the kallisto t2g file to your custom config file. 
### Create your own config file using the template:
The config file must be updated to the list where the top directory containing the pod5 folder and the pipeline outputs are located as well as the work directory. Ideally the config file would also be in the top directory folder. 
### Running dogme:
Running Dogme on typical dataset can take more than 24 hours, therefore it is recommended to run Dogme within a job or a saved virtual terminal such as screen or Tmux.  Change your folder to be where you want Dogme to run (the 'launch' directory), and launch Dogme directly from github using the following command: 
 ```
  nextflow run mortazavilab/dogme -c yourconfig.conf
  ```
By default, the pipeline will create several folders within the launch directory such as bams, bedMethyl, fastqs, and kallisto - all of which can be changed in the config file. If you need to resume your work add '-resume' to the nextflow command after deleting the html report and trace files.
