# dogme
A nextflow pipeline for basecalling nanopore reads with and without modifications before mapping and processing. Dogme supports 3 different read types: direct RNA ("RNA"), cDNA ("CDNA"), and genomic DNA ("DNA"). Given a folder of pod5 files, Dogme will download the latest dorado models, run dorado, merge bams, and map them with minimap2. For RNA and DNA, it will extract modifications using modkit and filter those. For RNA, cDNA, it will also extract a fastq and run kallisto.
### Installation:
The current version of Dogme has been tested on SLURM clusters and on Macs with the following software versions: 

| Software | Version |
|----------|---------|
| dorado   | 1.0     |
| samtools | 1.15.1  |
| minimap2 | 2.28    |
| mod_kit  | 0.4.3   |
| kallisto | 0.51.1  |
| bustools | 0.43.2  |

You must first install dorado, minimap2, samtools, modkit, kallisto, and bustools, and the latest version of nextflow and add their paths to the dogme.profile file in the launch directory or create an empty file of the same name. 
### Prerequisites:
Prerequisites include installing nextflow, using java/17 or better, installing modkit, installing minimap2, installing kallisto (compiled for long-reads), and installing bustools.
By default, nextflow will use the launchdirectory as the place to create its workfolder.
You must also create your own long-read kallisto index and t2g file for your genome of interest. You will need to add the path to the genome fasta, the transcriptome gtf, the kallisto index and the kallisto t2g file to your custom config file. 
### Create your own config file using the template:
The config file must be updated to the list where the top directory containing the pod5 folder and the pipeline outputs are located as well as the work directory. Ideally the config file would also be in the top directory folder. For most use cases only the top parameters in the config file will need to be changed. 

``` 
params {
    sample = 'C2C12_RNA_r1'
    //readType can either be 'RNA', 'DNA' or 'CDNA'
    readType = 'RNA'
    // modification type should be set as necessary if different from 'inosine_m6A,pseU,m5C' for RNA and '5mCG_5hmCG,6mA' for DNA. 
    //modifications = 'inosine_m6A,pseU,m5C'
    //change setting if necessary 
    minCov = 3
    perMod = 5
    // change if the launch directory is not where the pod5 and output directories should go
    topDir = "${launchDir}"

    // the following file should be edited to add all the necessary paths for commands such as
    // dorado, samtools, minimap2, kallisto, and bustools
    scriptEnv = "${launchDir}/dogme.profile"

    // needs to be modified to match the right genomic reference
     genome_annot_refs = [
     [name: 'mm39', genome: '/path/genomeRef/IGVFFI9282QLXO.fasta', annot: '/path/genomeRef/IGVFFI4777RDZK.gtf'],
     [name: 'C57BL_6J_T2T_v1', genome: '/path/genomeRef/C57BL_6J_T2T_v1/unmasked.fa', annot: '/path/genomeRef/C57BL_6J_T2T_v1/genes.gtf'],
     [name: 'CAST_EiJ_T2T_v1', genome: '/path/genomeRef/CAST_EiJ_T2T_v1/unmasked.fa', annot: '/path/genomeRef/CAST_EiJ_T2T_v1/genes.gtf']
     ]
    kallistoIndex = '/path/kallistoref/mm39GencM36_k63.idx'
    t2g = '/pathA/kallistoref/mm39GencM36_k63.t2g'
    
    //default accuracy is sup
    accuracy = "sup"
    // change this value if 0.9 is too strict
    // if set to null or '' then modkit will determine its threshold by sampling reads. 
    modkitFilterThreshold = 0.9

    //rest of config file - see dogmetest-param.conf
```
  
Be sure to change the process section of the example config file to reflect your cluster environment. 
### Running dogme:
Running Dogme on typical dataset can take more than 24 hours, therefore it is recommended to run Dogme within a job or a saved virtual terminal such as screen or Tmux.  Change your folder to be where you want Dogme to run (the 'launch' directory), and launch Dogme directly from github using the following command: 

 ```
  nextflow run mortazavilab/dogme -c yourconfig.conf
```
By default, the pipeline will create several folders within the launch directory such as bams, bedMethyl, fastqs, and kallisto - all of which can be customized in the config file. If you need to resume your work add '-resume' to the nextflow command after deleting the html report and trace files.

You can also use part of the workflow to only generate a report (-entry reports) or to remap the reads starting from unmapped bam file (-entry remap). 
