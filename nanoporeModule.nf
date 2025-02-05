#!/usr/bin/env nextflow


process softwareVTask {
    input:
    val version
    val modelPath
    output:
    path "${params.sample}.softwareVersion.txt"
    publishDir params.topDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    
    echo "dogme $version" > "${params.sample}.softwareVersion.txt"
    
    doradoV=\$(dorado -v 2>&1)
    echo "dorado \$doradoV" >> "${params.sample}.softwareVersion.txt"

    samtoolsV=\$(samtools version |grep samtools)
    echo \$samtoolsV >> "${params.sample}.softwareVersion.txt"

    minimap2V=\$(minimap2 --version 2>&1)
    echo "minimap2 \$minimap2V" >> "${params.sample}.softwareVersion.txt"

    modkitV=\$(modkit --version 2>&1)
    echo \$modkitV >> "${params.sample}.softwareVersion.txt"

    kallistoV=\$(kallisto version)
    echo \$kallistoV >> "${params.sample}.softwareVersion.txt"

    bustoolsV=\$(bustools version)
    echo \$bustoolsV >> "${params.sample}.softwareVersion.txt"
    
    echo "Dorado Models Used: " >> "${params.sample}.softwareVersion.txt"
    for folder in "${modelPath}"/*; do
    fullfile=\$(basename "\$folder")
    echo "\$fullfile" >> "${params.sample}.softwareVersion.txt"
    done
   
    """
}


process doradoDownloadTask {
    input:
    val dirPath
    val doradoModel
    output:
    val dirPath

    script:
    """
    echo "${dirPath}"
    . ${params.scriptEnv}
    mkdir -p ${dirPath}
    dorado download --data ${params.podDir} --model ${doradoModel}
    cp -rp *_* ${dirPath}
    """
}

process doradoTask {
    errorStrategy 'ignore'    
    input:
    path inputFile
    val modDirIgnore
    path modDirGood
    val doradoModel

    output:
    path "${inputFile.simpleName}.bam"
    publishDir params.dorDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    mkdir -p ${params.dorDir}

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

    publishDir params.bamDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    
     if [[ "${params.readType}" == "RNA" ]]; then
        minimap2_opts="-ax splice -uf --junc-bed ${params.annotRef}"
    elif [[ "${params.readType}" == "CDNA" ]]; then
        minimap2_opts="-ax splice:hq -uf --junc-bed ${params.annotRef}"  
    else
        minimap2_opts="-ax lr:hq"  
    fi
    
    
    samtools fastq --threads 64 -T MM,ML,pt ${params.sample}.unmapped.bam | \
    minimap2 -t 64 \$minimap2_opts --secondary=no --MD -y ${params.genomeRef} - | \
    samtools sort - --threads 64 > ${params.sample}.bam \
    && samtools index -@ 64 ${params.sample}.bam
    """
}

//splitbam files into plus and minus strands for direct rna
process separateStrandsTask {
    input: 
    path inputbam
    path inputbambai
    output:
    path "${params.sample}.plus.bam"
    path "${params.sample}.plus.bam.bai"
    path "${params.sample}.minus.bam"
    path "${params.sample}.minus.bam.bai"
    
    publishDir params.bamDir, mode: 'copy'
    
    script:
    """
    . ${params.scriptEnv}
    samtools view -b -f 16 ${inputbam} -o ${params.sample}.minus.bam && samtools index -@ 32 ${params.sample}.minus.bam
    samtools view -b -F 16 ${inputbam} -o ${params.sample}.plus.bam && samtools index -@ 32 ${params.sample}.plus.bam
    """
}

process modkitTask {
    input:
    path inputFile
    path inputFileBai
    
    output:
    path "${params.sample}.bed"
    
    publishDir params.bedDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    modkit pileup -t 12 --filter-threshold 0.9 ${inputFile} ${params.sample}.bed
    """
}

process modkitStrandedTask {
    input:
    path inputPlusFile
    path inputPlusFileBai
    path inputMinusFile
    path inputMinusFileBai
    
    output:
    path "${params.sample}.plus.bed"
    path "${params.sample}.minus.bed"
    
    publishDir params.bedDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    modkit pileup -t 12 --filter-threshold 0.9 ${inputPlusFile} ${params.sample}.plus.bed
    modkit pileup -t 12 --filter-threshold 0.9 ${inputMinusFile} ${params.sample}.minus.bed
    """
}

process filterbedTask {
    input:
    path inputFile

    output:
    path "${params.sample}.filtered-${params.minCov}-${params.perMod}.bed"
    
    publishDir params.bedDir, mode: 'copy'

    script:
    """
    python ${projectDir}/scripts/filterbed.py ${params.minCov} ${params.perMod} ${params.sample}.bed "${params.sample}.filtered-${params.minCov}-${params.perMod}.bed"
    """
}

process filterbedStrandedTask {
    input:
    path inputPlusFile
    path inputMinusFile
    output:
    path "${params.sample}.filtered-${params.minCov}-${params.perMod}.plus.bed"
    path "${params.sample}.filtered-${params.minCov}-${params.perMod}.minus.bed"
    
    publishDir params.bedDir, mode: 'copy'

    script:
    """
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

// The splitmodification task generates bed files for each DNA modification. 
// Modifications are identified by letters: 5mCG (m), 6mA (a),and hydroxymethylation (h). 
// The files generated by modkit are grepped for the letter codes. 
process splitModificationTask {
    input:
    path inputFile
    output:
    path "${params.sample}.5mCG.filtered.bed"
    path "${params.sample}.5hmCG.filtered.bed"
    path "${params.sample}.6mA.filtered.bed"
    publishDir params.bedDir, mode: 'copy'
    script:
    """
    # Extract 5mCG (methylation)
    grep -w 'm' ${inputFile} > "${params.sample}.5mCG.filtered.bed"
    # Extract 5hmCG (hydroxymethylation)
    grep -w 'h' ${inputFile} > "${params.sample}.5hmCG.filtered.bed"
    # Extract 6mA
    grep -w 'a' ${inputFile} > "${params.sample}.6mA.filtered.bed"
    """
}

// The splitmodification task generates bed files for each modification. 
// Modifications are identified by numbers and letters: inosine (17596), m5c (m), m6a (a), and pseU (17802). 
process splitModificationStrandedTask {
    input:
    path inputPlusFile
    path inputMinusFile
    output:
    path "${params.sample}.inosine.plus.filtered.bed"
    path "${params.sample}.inosine.minus.filtered.bed"
    path "${params.sample}.m6A.plus.filtered.bed"
    path "${params.sample}.m6A.minus.filtered.bed"
    path "${params.sample}.pseU.plus.filtered.bed"
    path "${params.sample}.pseU.minus.filtered.bed"
    path "${params.sample}.m5C.plus.filtered.bed"
    path "${params.sample}.m5C.minus.filtered.bed"
    publishDir params.bedDir, mode: 'copy'
    script:
    """
    # Extract m6A modifications (Plus & Minus strands)
    grep -w 'a' ${inputPlusFile} > "${params.sample}.m6A.plus.filtered.bed"
    grep -w 'a' ${inputMinusFile} > "${params.sample}.m6A.minus.filtered.bed"
    
    # Extract inosine modifications (Plus & Minus strands)
    grep -w '17596' ${inputPlusFile} > "${params.sample}.inosine.plus.filtered.bed"
    grep -w '17596' ${inputMinusFile} > "${params.sample}.inosine.minus.filtered.bed"

    # Extract pseudouridine (pseU) modifications (Plus & Minus strands)
    grep -w '17802' ${inputPlusFile} > "${params.sample}.pseU.plus.filtered.bed"
    grep -w '17802' ${inputMinusFile} > "${params.sample}.pseU.minus.filtered.bed"
    
    # Extract m5C modifications (Plus & Minus strands)
    grep -w 'm' ${inputPlusFile} > "${params.sample}.m5C.plus.filtered.bed"
    grep -w 'm' ${inputMinusFile} > "${params.sample}.m5C.minus.filtered.bed"
    """
}


workflow modWorkflow {
    take:
    theVersion
	theModel 
	modelDirectory
    
    main: 
	// Download the latest dorado models
     modelPath = doradoDownloadTask(modelDirectory, theModel)
    //Report all the software versions in report file  
    softwareVTask(theVersion, modelPath)
	def pod5FilesChannel = Channel.fromPath("${params.podDir}/*.pod5")
	// Run doradoTask for each input file
	bamFiles = doradoTask(pod5FilesChannel, modelPath, modelDirectory, theModel).collectFile()
	
	// Count all of the files as a way to force synchronization before merging
	fileCount = bamFiles.map { it.size() }.first()
	
	// Run merge task using the file count
	unmappedbam = mergeBamsTask(fileCount)
	
	// Run minimap
	mappedBams = minimapTask(unmappedbam)
	
    if (params.readType == 'RNA' || params.readType == 'CDNA') {
    // Run extractFastq
	fastqFile = extractfastqTask(unmappedbam)
	
	// Run kallistoTask using the extracted FASTQ file
	kallistoResults = kallistoTask(fastqFile)
    }
    
    if (params.readType == 'DNA') { 
    // DNA: Apply correct modifications (5mCG, 5hmCG, 6mA)
        bedfile = modkitTask(mappedBams)
        filterbed = filterbedTask(bedfile)
        splitResults = splitModificationTask(filterbed)
    }    
    
    if (params.readType == 'RNA') {
	// RNA: Apply stranded modifications (inosine_m6A, pseU, m5C)
        strandedBams = separateStrandsTask(mappedBams)
        bedfile = modkitStrandedTask(strandedBams)
        filterbed = filterbedStrandedTask(bedfile)
        splitResults = splitModificationStrandedTask(filterbed)
    }
}
