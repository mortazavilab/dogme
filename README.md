# dogme
A nextflow pipeline for modification calling in nanopore reads. 
### installation:
The current version of Dogme is customized to run on the UCI HPC3 cluster. In particular Dorado, Minimap, and Kallisto are already installed on hpc3. Future versions will be cluster independent.
### prerequisites:
Prerequisites include installing nextflow, using java/17 or better, installing dogme in an accessible folder, installing modkit, downloading dorado models into a defined folder. 
You must create your own work folder for nextflow on HPC3 outside of your home directory. 
### create your own config file using the template:
The config file must be updated to the list where the top directory containing the pod5 folder and the pipeline outputs are located as well as the work directory. Ideally the config file would also be in the top directory folder. 
### running dogme:
Running Dogme on typical dataset can take more than 24 hours, therefore it is recommended to run Dogme in Tmux mode with srun. After loading Java/17 run Dogme using the following command : 
  'nextflow run /path_to_dogme/dogme.nf -c dogmetest-param.conf'
  and if you need to resume your work add '-resume' to the nextflow command 
