# Be sure to include the path to dorado, minimap2, samtools, bedtools2, kallisto, and bustools in your PATH
#module load dorado
#module load minimap2
#module load samtools
#module load bedtools2
#module load kallisto
#export PATH=$PATH:/path/kb_python/kb_python/bins/linux/bustools
# the following should be set to point to where the latest version of modkit and the models are located
#export MODKITBASE=/share/crsp/lab/seyedam/share/modkit
#export MODKITMODEL=${MODKITBASE}/models/r1041_e82_400bps_hac_v5.2.0@v0.1.0
# the following is only needed if using the torch version of modkit
#export LIBTORCH=${MODKITBASE}/libtorch
#export DYLD_LIBRARY_PATH=${LIBTORCH}/lib
#export LD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}
#export PATH=${MODKITBASE}/dist_modkit_v0.5.0_5120ef7_tch:$PATH