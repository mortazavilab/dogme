#!/usr/bin/env nextflow

process doradoDownloadTask {
    input:
    val dirPath
    val doradoModel
    output:
    val "*_*/"
    publishDir { dirPath } , mode: 'symlink'

    script:
    """
    echo "${dirPath}"
    . ${params.scriptEnv}
    mkdir -p ${dirPath}
    dorado download --data ${params.podDir} --model ${doradoModel}
    """
}

process doradoTask {
    // errorStrategy 'ignore'    
    input:
    path inputFile
    val modDirIgnore
    val modDirGood
    val doradoModel

    output:
    path "${inputFile.simpleName}.bam"
    publishDir params.dorDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    mkdir -p ${params.dorDir}
    fullfile=\$(basename $inputFile)
    basefile=\${fullfile%.*}

    dorado basecaller ${doradoModel} --models-directory ${modDirGood}  --estimate-poly-a --batchsize 32 $inputFile > "${inputFile.simpleName}.bam"
    """
}

process mergeBamsTask {
    input:
    val fileCount
    output:
    path "${params.sample}.unmapped.bam"

    publishDir params.bamDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    samtools merge --threads ${task.cpus} -o ${params.sample}.unmapped.bam ${params.dorDir}/*.bam
    """
}

process minimapTask {
    input:
    path inputFile
    output:
    path "${params.sample}.bam"
    path "${params.sample}.bam.bai"
    path "${params.sample}.plus.bam"
    path "${params.sample}.plus.bam.bai"
    path "${params.sample}.minus.bam"
    path "${params.sample}.minus.bam.bai"

    publishDir params.bamDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    samtools fastq --threads 64 -T MM,ML,pt ${params.sample}.unmapped.bam | \
    minimap2 -t 64 -ax splice --junc-bed ${params.annotRef} --secondary=no --MD -y ${params.genomeRef} - | \
    samtools sort - --threads 64 > ${params.sample}.bam \
    && samtools index -@ 64 ${params.sample}.bam
    samtools view -b -f 16 ${params.sample}.bam -o ${params.sample}.minus.bam && samtools index -@ 32 ${params.sample}.minus.bam
    samtools view -b -F 16 ${params.sample}.bam -o ${params.sample}.plus.bam && samtools index -@ 32 ${params.sample}.plus.bam
    """
}

process modkitTask {
    input:
    path inputFile
    path inputFileBai
    path inputPlusFile
    path inputPlusFileBai
    path inputMinusFile
    path inputMinusFileBai
    
    output:
    path "${params.sample}.bed"
    path "${params.sample}.plus.bed"
    path "${params.sample}.minus.bed"    
    
    publishDir params.bedDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    modkit pileup -t 12 --filter-threshold 0.9 ${params.sample}.bam ${params.sample}.bed
    modkit pileup -t 12 --filter-threshold 0.9 ${params.sample}.plus.bam ${params.sample}.plus.bed
    modkit pileup -t 12 --filter-threshold 0.9 ${params.sample}.minus.bam ${params.sample}.minus.bed
    """
}

process filterbedTask {
    input:
    path inputFile
    path inputPlusFile
    path inputMinusFile
    output:
    path "${params.sample}.filtered-${params.minCov}-${params.perMod}.bed"
    path "${params.sample}.filtered-${params.minCov}-${params.perMod}.plus.bed"
    path "${params.sample}.filtered-${params.minCov}-${params.perMod}.minus.bed"
    
    publishDir params.bedDir, mode: 'copy'

    script:
    """
    python ${projectDir}/scripts/filterbed.py ${params.minCov} ${params.perMod} ${params.sample}.bed "${params.sample}.filtered-${params.minCov}-${params.perMod}.bed"
    python ${projectDir}/scripts/filterbed.py ${params.minCov} ${params.perMod} ${params.sample}.plus.bed "${params.sample}.filtered-${params.minCov}-${params.perMod}.plus.bed"
    python ${projectDir}/scripts/filterbed.py ${params.minCov} ${params.perMod} ${params.sample}.minus.bed "${params.sample}.filtered-${params.minCov}-${params.perMod}.minus.bed"
    """
}

process extractfastqTask {
    input:
    path inputFile
    output:
    path "${params.sample}.fastq.gz"

    publishDir params.fastqDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    samtools fastq --threads 6 ${params.sample}.unmapped.bam > ${params.sample}.fastq
    gzip -v ${params.sample}.fastq
    """
}

process kallistoTask {
    input:
    path inputFile
    output:
    path "${params.sample}"

    publishDir params.kallistoDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    mkdir -p ${params.sample}


    kallisto bus --long --threshold 0.8 -x bulk -i ${params.kallistoIndex} -t ${task.cpus} -o ${params.sample} "${inputFile}"

    bustools sort -t ${task.cpus} ${params.sample}/output.bus -o ${params.sample}/sorted.bus

    bustools count ${params.sample}/sorted.bus -t ${params.sample}/transcripts.txt  -e ${params.sample}/matrix.ec  -o ${params.sample}/count --cm -m -g ${params.t2g}

    kallisto quant-tcc -t ${task.cpus} --long -P ONT ${params.sample}/count.mtx -i ${params.kallistoIndex} -f ${params.sample}/flens.txt -e ${params.sample}/count.ec.txt -o ${params.sample}
    """
}

// The splitmodification task generates bed files for each modification. 
// Modifications are identified by numbers and letters: inosine (17596), m5c (m), m6a (a), and pseU (17802). 
// The files generated by modkit are grepped for the letter codes. 
// For the numerical codes, we first grepped for the code and then filtered out for the other codes to avoid 
// accidentally selecting records that match chromosome coordinates. 
process  splitModificationTask {
    input:
    path inputFile
    path inputPlusFile
    path inputMinusFile
    
    output:
    path "${params.sample}.m6a.filtered-${params.minCov}-${params.perMod}.bed"
    path "${params.sample}.inosine.filtered-${params.minCov}-${params.perMod}.bed"
    path "${params.sample}.m6a.filtered-${params.minCov}-${params.perMod}.plus.bed"
    path "${params.sample}.inosine.filtered-${params.minCov}-${params.perMod}.plus.bed"
    path "${params.sample}.m6a.filtered-${params.minCov}-${params.perMod}.minus.bed"
    path "${params.sample}.inosine.filtered-${params.minCov}-${params.perMod}.minus.bed"
    
    path "${params.sample}.m5c.filtered-${params.minCov}-${params.perMod}.bed"
    path "${params.sample}.pseU.filtered-${params.minCov}-${params.perMod}.bed"
    path "${params.sample}.m5c.filtered-${params.minCov}-${params.perMod}.plus.bed"
    path "${params.sample}.pseU.filtered-${params.minCov}-${params.perMod}.plus.bed"
    path "${params.sample}.m5c.filtered-${params.minCov}-${params.perMod}.minus.bed"
    path "${params.sample}.pseU.filtered-${params.minCov}-${params.perMod}.minus.bed"

    publishDir params.bedDir, mode: 'copy'

    script:
    """
    grep a ${inputFile} > "${params.sample}.m6a.filtered-${params.minCov}-${params.perMod}.bed"
    grep 17596 ${inputPlusFile} | grep -v a |grep -v m | grep -v 17802 ${inputFile} > "${params.sample}.inosine.filtered-${params.minCov}-${params.perMod}.bed"
    grep a ${inputPlusFile} > "${params.sample}.m6a.filtered-${params.minCov}-${params.perMod}.plus.bed"
    grep 17596 ${inputPlusFile} | grep -v a |grep -v m | grep -v 17802 ${inputFile} > "${params.sample}.inosine.filtered-${params.minCov}-${params.perMod}.plus.bed"
    grep a ${inputMinusFile} > "${params.sample}.m6a.filtered-${params.minCov}-${params.perMod}.minus.bed"
    grep 17596 ${inputPlusFile} | grep -v a |grep -v m | grep -v 17802 ${inputFile} > "${params.sample}.inosine.filtered-${params.minCov}-${params.perMod}.minus.bed"
    
    grep m ${inputFile} > "${params.sample}.m5c.filtered-${params.minCov}-${params.perMod}.bed"
    grep 17802 ${inputFile} | grep -v a |grep -v m | grep -v 17596 > "${params.sample}.pseU.filtered-${params.minCov}-${params.perMod}.bed"
    grep m ${inputPlusFile} > "${params.sample}.m5c.filtered-${params.minCov}-${params.perMod}.plus.bed"
    grep 17802 ${inputPlusFile} | grep -v a |grep -v m | grep -v 17596 > "${params.sample}.pseU.filtered-${params.minCov}-${params.perMod}.plus.bed"
    grep m ${inputMinusFile} > "${params.sample}.m5c.filtered-${params.minCov}-${params.perMod}.minus.bed"
    grep 17802 ${inputMinusFile} | grep -v a |grep -v m | grep -v 17596 > "${params.sample}.pseU.filtered-${params.minCov}-${params.perMod}.minus.bed"
    """
}


workflow modWorkflow {
    take:
	theModel 
	modelDirectory
    
    main: 
	// Download the latest dorado models
        modelPath = doradoDownloadTask(modelDirectory, theModel)
    
	def pod5FilesChannel = Channel.fromPath("${params.podDir}/*.pod5")
	// Run doradoTask for each input file
	bamFiles = doradoTask(pod5FilesChannel, modelPath, modelDirectory, theModel).collectFile()
	
	// Count all of the files as a way to force synchronization before merging
	fileCount = bamFiles.map { it.size() }.first()
	
	// Run merge task using the file count
	unmappedbam = mergeBamsTask(fileCount)
	
	// Run minimap
	mappedBams = minimapTask(unmappedbam)
	
	// Run extractFastq
	fastqFile = extractfastqTask(unmappedbam)
	
	// Run kallistoTask using the extracted FASTQ file
	kallistoResults = kallistoTask(fastqFile)
	
	// Run modkit
	bedfile = modkitTask(mappedBams)
	
	// Filter BED file
	filterbed = filterbedTask(bedfile) 
	
	// split the combined bed files into a bedfile for each modification 
	splitResults = splitModificationTask(filterbed)      
}




