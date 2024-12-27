#!/usr/bin/env nextflow

process doradoTask {
    // errorStrategy 'ignore'    
    input:
    path inputFile
    val rnaMod
    val rnaModel

    output:
    path "${inputFile.simpleName}.${rnaMod}.bam"
    publishDir params.dorDir, mode: 'copy'

    script:
    """
    . ${projectDir}/dogme.profile
    mkdir -p ${params.dorDir}
    fullfile=\$(basename $inputFile)
    basefile=\${fullfile%.*}

    dorado basecaller ${rnaModel} --models-directory ${params.rnaModelDir}  --estimate-poly-a --batchsize 32 $inputFile > "${inputFile.simpleName}.${rnaMod}.bam"
    """
}

process mergeBamsTask {
    input:
    val fileCount
    val rnaMod
    output:
    path "${params.sample}.${rnaMod}.unmapped.bam"

    publishDir params.bamDir, mode: 'copy'

    script:
    """
   . ${projectDir}/dogme.profile
    samtools merge --threads ${task.cpus} -o ${params.sample}.${rnaMod}.unmapped.bam ${params.dorDir}/*.${rnaMod}.bam
    chmod 775 "${params.sample}.${rnaMod}.unmapped.bam"
    """
}

process minimapTask {
    input:
    path inputFile
    val rnaMod
    output:
    path "${params.sample}.${rnaMod}.bam"
    path "${params.sample}.${rnaMod}.bam.bai"

    publishDir params.bamDir, mode: 'copy'

    script:
    """
    . ${projectDir}/dogme.profile
    samtools fastq --threads 64 -T MM,ML,pt ${params.sample}.${rnaMod}.unmapped.bam | \
    minimap2 -t 64 -ax splice --junc-bed ${params.annotRef} --secondary=no --MD -y ${params.genomeRef} - | \
    samtools sort - --threads 64 > ${params.sample}.${rnaMod}.bam \
    && samtools index -@ 64 ${params.sample}.${rnaMod}.bam
    chmod 775 "${params.sample}.${rnaMod}.bam"
    chmod 775 "${params.sample}.${rnaMod}.bam.bai"
    """
}

process modkitTask {
    input:
    path inputFile
    path inputFileBai
    val rnaMod
    output:
    path "${params.sample}.${rnaMod}.bed"

    publishDir params.bedDir, mode: 'copy'

    script:
    """
    . ${projectDir}/dogme.profile
    modkit pileup -t 12 --filter-threshold 0.9 ${params.sample}.${rnaMod}.bam ${params.sample}.${rnaMod}.bed
    chmod 775 "${params.sample}.${rnaMod}.bed"
    """
}

process filterbedTask {
    input:
    path inputFile
    val rnaMod
    output:
    path "${params.sample}.${rnaMod}.filtered-${params.minCov}-${params.perMod}.bed"

    publishDir params.bedDir, mode: 'copy'

    script:
    """
    python ${projectDir}/scripts/filterbed.py ${params.minCov} ${params.perMod} ${params.sample}.${rnaMod}.bed "${params.sample}.${rnaMod}.filtered-${params.minCov}-${params.perMod}.bed"
    chmod 775 "${params.sample}.${rnaMod}.filtered-${params.minCov}-${params.perMod}.bed"
    """
}

process extractfastqTask {
    input:
    path inputFile
    val rnaMod
    output:
    path "${params.sample}.${rnaMod}.fastq.gz"

    publishDir params.fastqDir, mode: 'copy'

    script:
    """
    . ${projectDir}/dogme.profile
    samtools fastq --threads 6 ${params.sample}.${rnaMod}.unmapped.bam > ${params.sample}.${rnaMod}.fastq
    gzip -v ${params.sample}.${rnaMod}.fastq
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
    . ${projectDir}/dogme.profile
    mkdir -p ${params.sample}


    kallisto bus --long --threshold 0.8 -x bulk -i ${params.kallistoIndex} -t ${task.cpus} -o ${params.sample} "${inputFile}"

    bustools sort -t ${task.cpus} ${params.sample}/output.bus -o ${params.sample}/sorted.bus

    bustools count ${params.sample}/sorted.bus -t ${params.sample}/transcripts.txt  -e ${params.sample}/matrix.ec  -o ${params.sample}/count --cm -m -g ${params.t2g}

    kallisto quant-tcc -t ${task.cpus} --long -P ONT ${params.sample}/count.mtx -i ${params.kallistoIndex} -f ${params.sample}/flens.txt -e ${params.sample}/count.ec.txt -o ${params.sample}

    chmod 775 "${params.sample}"
    """
}


process  splitM6aInosineTask {
    input:
    path inputFile
    output:
    path "${params.sample}.m6a.filtered-${params.minCov}-${params.perMod}.bed"
    path "${params.sample}.inosine.filtered-${params.minCov}-${params.perMod}.bed"

    publishDir params.bedDir, mode: 'copy'

    script:
    """
    grep a ${inputFile} > "${params.sample}.m6a.filtered-${params.minCov}-${params.perMod}.bed"
    grep -v a ${inputFile} > "${params.sample}.inosine.filtered-${params.minCov}-${params.perMod}.bed"
    chmod 775 "${params.sample}.m6a.filtered-${params.minCov}-${params.perMod}.bed" "${params.sample}.inosine.filtered-${params.minCov}-${params.perMod}.bed"
    """
}


workflow modWorkflow {
    take:
	theMod
	theModel 
    main: 
	println "theMod: " + theMod 
	def pod5FilesChannel = Channel.fromPath("${params.podDir}/*.pod5")
    // Run doradoTask for each input file
    bamFiles = doradoTask(pod5FilesChannel, theMod, theModel).collectFile()

    // Count all of the files as a way to force synchronization before merging
    fileCount = bamFiles.map { it.size() }.first()

    // Run merge task using the file count
    unmappedbam = mergeBamsTask(fileCount, theMod)

    // Run minimap
    mappedBams = minimapTask(unmappedbam, theMod)
    
    // Run extractFastq
    if (theMod == 'm5c') {
        fastqFile = extractfastqTask(unmappedbam, theMod)

        // Run kallistoTask using the extracted FASTQ file
        kallistoResults = kallistoTask(fastqFile)
    }

    // Run modkit
    bedfile = modkitTask(mappedBams, theMod)

    // Filter BED file
    filterbed = filterbedTask(bedfile, theMod)
    
    if (theMod == 'm6aInosine') {
        splitResults = splitM6aInosineTask(filterbed)
        }
}

// Export the workflow
//export modWorkflow
