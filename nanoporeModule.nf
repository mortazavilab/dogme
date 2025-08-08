#!/usr/bin/env nextflow
nextflow.enable.dsl=2


// Use Nextflow's 'include' to alias openChromatinTask for bg and bed
include { openChromatinTask as openChromatinTaskBg } from './modkitModule'
include { openChromatinTask as openChromatinTaskBed } from './modkitModule'

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
    
    if [[ "${params.readType}" == "DNA" ]] || [[ "${params.readType}" == "RNA" ]]; then
    modkitV=\$(modkit --version 2>&1)
    echo \$modkitV >> "${params.sample}.softwareVersion.txt"
    fi
    
    if [[ "${params.readType}" == "CDNA" ]] || [[ "${params.readType}" == "RNA" ]]; then
    kallistoV=\$(kallisto version)
    echo \$kallistoV >> "${params.sample}.softwareVersion.txt"
    bustoolsV=\$(bustools version)
    echo \$bustoolsV >> "${params.sample}.softwareVersion.txt"
    fi
    
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
    tuple path(inputFile), val(genomeRef), val(annotRef), val(genomeName)
    output:
    tuple path("${params.sample}.${genomeName}.bam"), path("${params.sample}.${genomeName}.bam.bai"), val(genomeName)
    publishDir params.bamDir, mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    if [[ "${params.readType}" == "RNA" ]]; then
        minimap2_opts="-ax splice -uf --junc-bed ${annotRef}"
    elif [[ "${params.readType}" == "CDNA" ]]; then
        minimap2_opts="-ax splice:hq -uf --junc-bed ${annotRef}"  
    else
        minimap2_opts="-ax lr:hq"  
    fi
    samtools fastq --threads 64 -T MM,ML,pt ${params.sample}.unmapped.bam | \
    minimap2 -t 64 \$minimap2_opts -L --secondary=no --MD -y ${genomeRef} - | \
    samtools sort - --threads 64 > ${params.sample}.${genomeName}.bam \
    && samtools index -@ 64 ${params.sample}.${genomeName}.bam
    """
}

process getChrListTask {
    tag "${genomeName}"

    input:
    tuple path(mapped_bam), path(mapped_bai), val(genomeName)

    output:
    tuple path(chr_sizes_file), path(mapped_bam), path(mapped_bai), val(genomeName)

    publishDir path: "${params.topDir}/stats", pattern: "${chr_sizes_file}", mode: 'copy'

    script:
    chr_sizes_file = "${params.sample}.${genomeName}.chr_sizes.txt"
    """
     . ${params.scriptEnv}   
    # Use cut -f 1,2 to get chromosome name and size (length)
    # Filter out unmapped reads (*) and alternative contigs
    samtools idxstats ${mapped_bam} | cut -f 1,2 | grep -v '*' > "${chr_sizes_file}"
    """
}


//splitbam files into plus and minus strands for direct rna
process separateStrandsTask {
    input:
    tuple path(inputbam), val(genomeName)
    output:
    tuple path("${params.sample}.${genomeName}.plus.bam"), path("${params.sample}.${genomeName}.plus.bam.bai"), val(genomeName), emit: plus_strand
    tuple path("${params.sample}.${genomeName}.minus.bam"), path("${params.sample}.${genomeName}.minus.bam.bai"), val(genomeName), emit: minus_strand

    publishDir params.bamDir, mode: 'copy'

    script:
    """
    . ${params.scriptEnv}
    samtools view -b -f 16 ${inputbam} -o ${params.sample}.${genomeName}.minus.bam && samtools index -@ 32 ${params.sample}.${genomeName}.minus.bam
    samtools view -b -F 16 ${inputbam} -o ${params.sample}.${genomeName}.plus.bam && samtools index -@ 32 ${params.sample}.${genomeName}.plus.bam
    """
}

process modkitTask {
    input:
    tuple path(inputFile), path(inputBai), val(genomeName)

    output:
    path "*.bed"

    publishDir params.bedDir, mode: 'copy'
    script:
    // Build the filter threshold argument conditionally
    def filterThresholdArg = ''
    if (params.modkitFilterThreshold != null && params.modkitFilterThreshold != '') {
        filterThresholdArg = "--filter-threshold ${params.modkitFilterThreshold}"
    }
    """
    . ${params.scriptEnv}
    bedFileOutput="${params.sample}.${genomeName}.bed" # Default value
    if [[ "${params.readType}" == "RNA" ]]; then
        if [[ "${inputFile}" == *".plus."* ]]; then
            bedFileOutput="${params.sample}.${genomeName}.plus.bed"
        elif [[ "${inputFile}" == *".minus."* ]]; then
            bedFileOutput="${params.sample}.${genomeName}.minus.bed"
        fi
    fi
    modkit pileup -t 12 ${filterThresholdArg} "${inputFile}" "\${bedFileOutput}"
    """
}

process filterbedTask {
    input:
    path inputFile
    output:
    path "*.filtered*.bed"
    publishDir params.bedDir, mode: 'copy'
    script:
    """
    output_prefix="${inputFile.baseName}"
    bedFileOutput="\${output_prefix}.filtered-${params.minCov}-${params.perMod}.bed"
    python ${projectDir}/scripts/filterbed.py ${params.minCov} ${params.perMod} "${inputFile}" \${bedFileOutput}
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
    path "*.filtered.*"
    publishDir params.bedDir, mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    if [[ "${params.readType}" == "DNA" ]]; then 
        # Extract 5mCG (methylation)
        if grep -q 'm' "${inputFile}"; then
            grep -w 'm' "${inputFile}" > "${inputFile.baseName}.5mCG.filtered.bed"
        fi
        # Extract 5hmCG (hydroxymethylation)
        if grep -q 'h' "${inputFile}"; then
            grep -w 'h' "${inputFile}" > "${inputFile.baseName}.5hmCG.filtered.bed"
        fi
        # Extract 6mA
        if grep -q 'a' "${inputFile}"; then
            grep -w 'a' "${inputFile}" > "${inputFile.baseName}.6mA.filtered.bed"
        fi
    elif [[ "${params.readType}" == "RNA" ]]; then
        base_name="\$(basename "${inputFile}" .bed)"
        # Extract m6A modifications (Plus & Minus strands)
        if grep -q 'a' "${inputFile}"; then
            grep -w 'a' "${inputFile}" > "\${base_name/filtered*/m6A.filtered}.bed"
        fi
        # Extract inosine modifications (Plus & Minus strands)
        if grep -q '17596' "${inputFile}"; then
            grep -w '17596' "${inputFile}" > "\${base_name/filtered*/inosine.filtered}.bed"
        fi
        # Extract pseudouridine (pseU) modifications (Plus & Minus strands)
        if grep -q '17802' "${inputFile}"; then
            grep -w '17802' "${inputFile}" > "\${base_name/filtered*/pseU.filtered}.bed"
        fi
        # Extract m5C modifications (Plus & Minus strands)
        if grep -q 'm' "${inputFile}"; then
            grep -w 'm' "${inputFile}" > "\${base_name/filtered*/m5C.filtered}.bed"
        fi
        # Extract Nm modifications (Plus & Minus strands)
        if grep -qE '19228|19229|19227|69426' "${inputFile}"; then
            grep -Ew '19228|19229|19227|69426' "${inputFile}" > "\${base_name/filtered*/Nm.filtered}.bed"
        fi
    fi
    """
}

process generateReport {
    tag "Generate metadata report"
    input:
    path report_inputs
    path results
    path openBeds
    output:
    path "report.tsv", emit: report
    publishDir params.topDir, mode: 'copy'
    script:
    """
    python ${projectDir}/scripts/generate_report.py -i ${report_inputs} -o report.tsv
    """
}

process consolidateOpenChromatinBedTask {
    tag "${params.sample}.${genomeName}"
    input:
    tuple val(genomeName), path(bed_files)
    output:
    path "${params.sample}.${genomeName}.m6Aopen.bed"
    publishDir "${params.topDir}/openChromatin", mode: 'copy'
    script:
    """
    cat \$(ls ${bed_files} | sort) > ${params.sample}.${genomeName}.m6Aopen.bed
    """
}

process consolidateOpenChromatinBgTask {
    tag "${params.sample}.${genomeName}"
    input:
    tuple val(genomeName), path(bg_files)
    output:
    path "${params.sample}.${genomeName}.m6Aopen.bg"
    publishDir "${params.topDir}/openChromatin", mode: 'copy'
    script:
    """
    cat \$(ls ${bg_files} | sort) | sort -k1,1 -k2,2n > ${params.sample}.${genomeName}.m6Aopen.bg
    """
}


workflow modificationWorkflow {
    take:
    mapped_bams_ch // Channel containing tuples: [path(bam), path(bai), val(genomeName)]
    model_name     // The 'theModel' variable, needed for the open chromatin check

    main:
    if (params.readType == 'RNA') {
        mappedBamsTuplesRNA = mapped_bams_ch.map { it -> tuple(*it) }
        mappedBamsForStrands = mappedBamsTuplesRNA.map { bam, bai, genomeName -> tuple(bam, genomeName) }
        strands = separateStrandsTask(mappedBamsForStrands)
        plusStrand = strands.plus_strand
        minusStrand = strands.minus_strand
        combinedStrand = plusStrand.concat(minusStrand)
        bedfiles = modkitTask(combinedStrand)
    } else if (params.readType == 'DNA') {
        mappedBamsTuplesDNA = mapped_bams_ch.map { it -> tuple(*it) }

        if (model_name.contains('6mA')) {
            chrListOutput = getChrListTask(mappedBamsTuplesDNA)
            chrInfoForOpenChromatinBed = chrListOutput.flatMap { chr_sizes_file, bam_file, bai_file, genomeName ->
                chr_sizes_file.readLines().collect { line ->
                    def parts = line.split('\\t')
                    tuple(parts[0], parts[1], bam_file, bai_file, genomeName, 'bed')
                }
            }
            chrInfoForOpenChromatinBg = chrListOutput.flatMap { chr_sizes_file, bam_file, bai_file, genomeName ->
                chr_sizes_file.readLines().collect { line ->
                    def parts = line.split('\\t')
                    tuple(parts[0], parts[1], bam_file, bai_file, genomeName, 'bg')
                }
            }
            openChromatinResultsBg = openChromatinTaskBg(chrInfoForOpenChromatinBg)
            openChromatinResultsBed = openChromatinTaskBed(chrInfoForOpenChromatinBed)

            // Collect all per-chromosome bed files per genome and consolidate
            openChromatinResultsBed
                .map { bedfile -> 
                    def genomeName = bedfile.getParent().getParent().getName()
                    tuple(genomeName, bedfile)
                }
                .groupTuple()
                .set { perGenomeBedFiles }

            consolidatedBeds = consolidateOpenChromatinBedTask(perGenomeBedFiles)

            // Collect all per-chromosome bg files per genome and consolidate
            openChromatinResultsBg
                .map { bgfile ->
                    def genomeName = bgfile.getParent().getParent().getName()
                    tuple(genomeName, bgfile)
                }
                .groupTuple()
                .set { perGenomeBgFiles }

            consolidatedBgs = consolidateOpenChromatinBgTask(perGenomeBgFiles)
        }

        mappedBamsForModkit = mappedBamsTuplesDNA.map { bam, bai, genomeName -> tuple(bam, bai, genomeName) }
        bedfiles = modkitTask(mappedBamsForModkit)
    }

    filterbeds = filterbedTask(bedfiles)
    splitResults = splitModificationTask(filterbeds)
    splitResultsReport = splitResults.last()
    
    if (model_name.contains('6mA')) {
        consolidatedBedReport = consolidatedBeds.last()
    } else {
        consolidatedBedReport = Channel.empty()
    }

    generateReport(launchDir, splitResultsReport, consolidatedBedReport)
}

workflow mainWorkflow {
    take:
    theVersion
    theModel 
    modelDirectory
    
    main: 
    modelPath = doradoDownloadTask(modelDirectory, theModel)
    softwareVTask(theVersion, modelPath)
    def pod5FilesChannel = Channel.fromPath("${params.podDir}/*.pod5")
    bamFiles = doradoTask(pod5FilesChannel, modelPath, modelDirectory, theModel).collectFile()
    fileCount = bamFiles.map { it.size() }.first()
    unmappedbam = mergeBamsTask(fileCount)
    
    if (params.readType == 'RNA' || params.readType == 'CDNA') {
        // Run extractFastq
        fastqFile = extractfastqTask(unmappedbam)
        // Run kallistoTask using the extracted FASTQ file
        kallistoResults = kallistoTask(fastqFile)
    }

    def genomeAnnotChannel = Channel.fromList(params.genome_annot_refs)
    unmappedBams = unmappedbam.combine(genomeAnnotChannel).map { bam, ref ->
        tuple(bam, ref.genome, ref.annot, ref.name)
    }
    mappedBams = minimapTask(unmappedBams)

    if (params.readType == 'RNA' || params.readType == 'DNA') {
        modificationWorkflow(mappedBams, theModel)
    } else {
        placeholder1 = Channel.empty()
        placeholder2 = Channel.empty()
        generateReport(launchDir, placeholder1, placeholder2)
    }
}

workflow remapWorkflow {
    take:
    theVersion
    theModel 
    modelDirectory

    main:
    def unmappedbam = Channel.fromPath("${params.bamDir}/*.unmapped.bam")
    def genomeAnnotChannel = Channel.fromList(params.genome_annot_refs)
    unmappedBams = unmappedbam.combine(genomeAnnotChannel).map { bam, ref ->
        tuple(bam, ref.genome, ref.annot, ref.name)
    }
    mappedBams = minimapTask(unmappedBams)

    if (params.readType == 'RNA' || params.readType == 'DNA') {
        modificationWorkflow(mappedBams, theModel)
    } else {
        placeholder1 = Channel.empty()
        placeholder2 = Channel.empty()
        generateReport(launchDir, placeholder1, placeholder2)
    }
}

workflow reportsWorkflow {
    take:
    theVersion
    modelDirectory
    
    main:
    softwareVTask(theVersion, modelDirectory)
    
    placeholder1 = Channel.empty()
    placeholder2 = Channel.empty()
    generateReport(launchDir, placeholder1, placeholder2)
}
