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

    python ${projectDir}/scripts/software_versions.py \\
        --version "${version}" \\
        --read-type "${params.readType}" \\
        --model-path "${modelPath}" \\
        --output "${params.sample}.softwareVersion.txt" \\
        --sample "${params.sample}"
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
    cp -r *_* ${dirPath}
    """
}

process doradoTask {
    errorStrategy { task.attempt < 4 ? 'retry' : 'ignore' }    
    input:
    path inputFile
    val modDirIgnore
    path modDirGood
    val doradoModel
    output:
    path "${inputFile.simpleName}.bam"
    publishDir params.dorDir, mode: 'copy'
    script:
    def kitNameArg = params.kitName ? "--kit-name ${params.kitName}" : ''
    """
    . ${params.scriptEnv}
    mkdir -p ${params.dorDir}
    dorado basecaller ${doradoModel} --models-directory ${modDirGood} ${kitNameArg} --estimate-poly-a --batchsize 32 $inputFile > "${inputFile.simpleName}.bam"
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

process doradoDemuxTask {
    tag "${params.sample} demultiplex"
    input:
    path inputBam
    output:
    path "demux/**/*.bam", emit: demuxed_bams
    publishDir "${params.dorDir}/demux", mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    mkdir -p demux
    dorado demux --no-classify --output-dir demux ${inputBam}
    """
}

process demuxSummaryTask {
    tag "${params.sample} demux summary"
    input:
    path inputBams
    output:
    path "${params.sample}.demux_summary.tsv"
    publishDir params.topDir, mode: 'copy'
    script:
    def keepCount = params.keepBarcodes ?: 0
    """
    . ${params.scriptEnv}
    : > raw_counts.tsv
    while IFS= read -r -d '' inputBam; do
        inputName="\${inputBam##*/}"
        inputName="\$(printf '%s' "\${inputName}" | tr '[:upper:]' '[:lower:]')"
        if [[ "\${inputName}" =~ barcode([0-9]+) ]]; then
            category="bc\${BASH_REMATCH[1]}"
        elif [[ "\${inputName}" == *unclassified* ]]; then
            category="unclassified"
        elif [[ "\${inputName}" == *no_barcode* || "\${inputName}" == *nobarcode* ]]; then
            category="no_barcode"
        else
            category="other"
        fi
        printf '%s\\t%s\\n' "\${category}" "\$(samtools view -c "\${inputBam}")" >> raw_counts.tsv
    done < <(find -L . -type f -name '*.bam' -print0)

    awk -F '\\t' '{ counts[\$1] += \$2 } END { for (category in counts) print category "\\t" counts[category] }' raw_counts.tsv \
        | sort -k2,2nr -k1,1 > all_counts.tsv

    awk '\$1 ~ /^bc[0-9]+\$/ { print }' all_counts.tsv \
        | head -n ${keepCount} \
        | cut -f1 > selected_barcodes.txt
    if [[ ${keepCount} -eq 0 ]]; then
        awk '\$1 ~ /^bc[0-9]+\$/ { print \$1 }' all_counts.tsv > selected_barcodes.txt
    fi

    totalReads="\$(awk -F '\\t' '{ total += \$2 } END { print total + 0 }' all_counts.tsv)"
    awk -F '\\t' -v totalReads="\${totalReads}" '
        FNR == NR { selected[\$1] = 1; next }
        {
            retained = selected[\$1] ? "yes" : "no"
            percent = totalReads ? (100 * \$2 / totalReads) : 0
            printf "%s\\t%s\\t%.2f\\t%s\\n", \$1, \$2, percent, retained
        }
    ' selected_barcodes.txt all_counts.tsv > summary_rows.tsv

    {
        printf 'category\\treads\\tpercent_of_demuxed_reads\\tretained_for_downstream\\n'
        cat summary_rows.tsv
    } > ${params.sample}.demux_summary.tsv
    """
}

process selectTopDemuxBamsTask {
    tag "top ${params.keepBarcodes} barcodes"
    input:
    path inputBams
    output:
    path "selected/*.bam", emit: selected_bams
    script:
    """
    . ${params.scriptEnv}
    : > read_counts.tsv
    while IFS= read -r -d '' inputBam; do
        if [[ "\${inputBam}" =~ barcode([0-9]+) ]]; then
            barcode="\${BASH_REMATCH[1]}"
            readCount="\$(samtools view -c "\${inputBam}")"
            printf '%s\\t%s\\t%s\\n' "\${barcode}" "\${readCount}" "\${inputBam}" >> read_counts.tsv
        fi
    done < <(find -L . -type f -name '*.bam' -print0)

    awk -F '\\t' '{ counts[\$1] += \$2 } END { for (barcode in counts) print barcode "\\t" counts[barcode] }' read_counts.tsv \
        | sort -k2,2nr -k1,1n \
        | head -n ${params.keepBarcodes} \
        | cut -f1 > selected_barcodes.txt

    mkdir selected
    while IFS= read -r barcode; do
        mapfile -t barcodeBams < <(awk -F '\\t' -v barcode="\${barcode}" '\$1 == barcode { print \$3 }' read_counts.tsv)
        if [[ "\${#barcodeBams[@]}" -eq 1 ]]; then
            cp "\${barcodeBams[0]}" "selected/bc\${barcode}.bam"
        else
            samtools merge --threads ${task.cpus} -o "selected/bc\${barcode}.bam" "\${barcodeBams[@]}"
        fi
    done < selected_barcodes.txt
    """
}

process prepareDemuxBamTask {
    tag "${params.sample}.${barcode}"
    input:
    tuple path(inputBam), val(barcode)
    output:
    path "${params.sample}.${barcode}.unmapped.bam", emit: classified_unmapped
    publishDir params.bamDir, mode: 'copy'
    script:
    """
    cp ${inputBam} ${params.sample}.${barcode}.unmapped.bam
    """
}

process samtoolsImportTask {
    input:
    path inputFastq
    output:
    path "${params.sample}.unmapped.bam"
    publishDir params.bamDir, mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    samtools import -0 ${inputFastq} -o ${params.sample}.unmapped.bam
    """
}

process minimapTask {
    input:
    tuple path(inputFile), val(genomeRef), val(annotRef), val(genomeName)
    output:
    tuple path("${inputFile.simpleName}.${genomeName}.bam"), path("${inputFile.simpleName}.${genomeName}.bam.bai"), val(genomeName), emit: mapped_bams
    publishDir params.bamDir, mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    # Extract sample name from inputFile
    # Use the inputFile name to determine the sample name
    sample_name=${inputFile.simpleName}
    
    if [[ "${params.readType}" == "RNA" ]]; then
        python ${projectDir}/scripts/gtf_to_junction_bed.py ${annotRef} > junc.bed
        minimap2_opts="-ax splice -uf -G 500000 --junc-bed junc.bed"
    elif [[ "${params.readType}" == "CDNA" ]]; then
        python ${projectDir}/scripts/gtf_to_junction_bed.py ${annotRef} > junc.bed
        minimap2_opts="-ax splice:hq -uf -G 500000 --junc-bed junc.bed"
    else
        minimap2_opts="-ax lr:hq"
    fi
    samtools fastq --threads 64 -T MM,ML,pt \${sample_name}.unmapped.bam > \${sample_name}.fastq
    minimap2 -t 64 \$minimap2_opts -L --secondary=no --MD -y ${genomeRef} \${sample_name}.fastq > \${sample_name}.${genomeName}.sam
    rm \${sample_name}.fastq 
    samtools sort \${sample_name}.${genomeName}.sam --threads 64 > \${sample_name}.${genomeName}.bam \
    && samtools index -@ 64 \${sample_name}.${genomeName}.bam
    rm \${sample_name}.${genomeName}.sam
    """
}

process demuxMinimapTask {
    tag "${sampleName}.${genomeName}"
    input:
    tuple path(inputFile), val(sampleName), val(genomeRef), val(annotRef), val(genomeName)
    output:
    tuple path("${sampleName}.${genomeName}.bam"), path("${sampleName}.${genomeName}.bam.bai"), val(genomeName), val(sampleName), emit: mapped_bams
    publishDir params.bamDir, mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    if [[ "${params.readType}" == "RNA" ]]; then
        python ${projectDir}/scripts/gtf_to_junction_bed.py ${annotRef} > junc.bed
        minimap2_opts="-ax splice -uf -G 500000 --junc-bed junc.bed"
    elif [[ "${params.readType}" == "CDNA" ]]; then
        python ${projectDir}/scripts/gtf_to_junction_bed.py ${annotRef} > junc.bed
        minimap2_opts="-ax splice:hq -uf -G 500000 --junc-bed junc.bed"
    else
        minimap2_opts="-ax lr:hq"
    fi
    samtools fastq --threads 64 -T MM,ML,pt ${inputFile} > ${sampleName}.fastq
    minimap2 -t 64 \$minimap2_opts -L --secondary=no --MD -y ${genomeRef} ${sampleName}.fastq > ${sampleName}.${genomeName}.sam
    rm ${sampleName}.fastq
    samtools sort ${sampleName}.${genomeName}.sam --threads 64 > ${sampleName}.${genomeName}.bam \
    && samtools index -@ 64 ${sampleName}.${genomeName}.bam
    rm ${sampleName}.${genomeName}.sam
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

process demuxSeparateStrandsTask {
    input:
    tuple path(inputbam), val(genomeName), val(sampleName)
    output:
    tuple path("${sampleName}.${genomeName}.plus.bam"), path("${sampleName}.${genomeName}.plus.bam.bai"), val(genomeName), val(sampleName), emit: plus_strand
    tuple path("${sampleName}.${genomeName}.minus.bam"), path("${sampleName}.${genomeName}.minus.bam.bai"), val(genomeName), val(sampleName), emit: minus_strand
    publishDir params.bamDir, mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    samtools view -b -f 16 ${inputbam} -o ${sampleName}.${genomeName}.minus.bam && samtools index -@ 32 ${sampleName}.${genomeName}.minus.bam
    samtools view -b -F 16 ${inputbam} -o ${sampleName}.${genomeName}.plus.bam && samtools index -@ 32 ${sampleName}.${genomeName}.plus.bam
    """
}

process modkitTask {
    input:
    tuple path(inputFile), path(inputBai), val(genomeName)

    output:
    path "*.bed.gz"

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
    gzip "\${bedFileOutput}"
    """
}

process demuxModkitTask {
    input:
    tuple path(inputFile), path(inputBai), val(genomeName), val(sampleName)
    output:
    path "*.bed.gz", optional: true
    publishDir params.bedDir, mode: 'copy'
    script:
    def filterThresholdArg = ''
    if (params.modkitFilterThreshold != null && params.modkitFilterThreshold != '') {
        filterThresholdArg = "--filter-threshold ${params.modkitFilterThreshold}"
    }
    """
    . ${params.scriptEnv}
    bedFileOutput="${sampleName}.${genomeName}.bed"
    if [[ "${params.readType}" == "RNA" ]]; then
        if [[ "${inputFile}" == *".plus."* ]]; then
            bedFileOutput="${sampleName}.${genomeName}.plus.bed"
        elif [[ "${inputFile}" == *".minus."* ]]; then
            bedFileOutput="${sampleName}.${genomeName}.minus.bed"
        fi
    fi
    if [[ "\$(samtools view -c -F 4 ${inputFile})" -gt 0 ]]; then
        modkit pileup -t 12 ${filterThresholdArg} "${inputFile}" "\${bedFileOutput}"
        gzip "\${bedFileOutput}"
    fi
    """
}

process filterbedTask {
    input:
    path inputFile
    output:
    path "*.filtered*.bed.gz"
    publishDir params.bedDir, mode: 'copy'
    script:
    """
    output_prefix="\$(basename "${inputFile}" .bed.gz)"
    bedFileOutput="\${output_prefix}.filtered-${params.minCov}-${params.perMod}.bed.gz"
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

process demuxExtractfastqTask {
    tag "${sampleName} FASTQ"
    input:
    tuple path(inputFile), val(sampleName)
    output:
    tuple path("${sampleName}.fastq.gz"), val(sampleName)
    publishDir params.fastqDir, mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    samtools fastq --threads 6 ${inputFile} > ${sampleName}.fastq
    gzip -v ${sampleName}.fastq
    """
}

process demuxGenerateSeqspecTask {
    tag "${sampleName} seqspec"
    container 'ghcr.io/mortazavilab/dogme-pipeline:latest'
    input:
    tuple path(inputFastq), val(sampleName)
    path templates
    path scripts
    output:
    tuple path("${sampleName}.seqspec.yaml"), val(sampleName)
    publishDir params.fastqDir, mode: 'copy'
    script:
    def outputSpec = "${sampleName}.seqspec.yaml"
    def singleCellEnabled = params.containsKey('singleCell') && params.singleCell
    def singleCellKitArg = params.containsKey('singleCellKit') && params.singleCellKit ? "--single-cell-kit ${params.singleCellKit}" : ''
    def seqspecMd5Enabled = !params.containsKey('seqspecMd5') || params.seqspecMd5
    def seqspecTemplateArg = params.containsKey('seqspecTemplate') && params.seqspecTemplate ? "--template ${params.seqspecTemplate}" : ''
    def seqspecVariablesArg = params.containsKey('seqspecVariables') && params.seqspecVariables ? "--variables ${params.seqspecVariables}" : ''
    """
    . ${params.scriptEnv}
    python ${scripts}/generate_seqspec.py \
        --template-dir ${templates}/seqspec \
        --read-type ${params.readType} \
        ${singleCellEnabled ? '--single-cell' : ''} \
        ${singleCellKitArg} \
        ${seqspecMd5Enabled ? '' : '--no-md5'} \
        ${seqspecTemplateArg} \
        ${seqspecVariablesArg} \
        --fastq ${inputFastq} \
        --output ${outputSpec}
    """
}

process generateSeqspecTask {
    tag "${params.sample} seqspec"
    container 'ghcr.io/mortazavilab/dogme-pipeline:latest'
    input:
    path inputFastq
    path templates
    path scripts
    output:
    path "*.seqspec.yaml"
    publishDir params.fastqDir, mode: 'copy'
    script:
    def outputSpec = "${params.sample}.seqspec.yaml"
    def singleCellEnabled = params.containsKey('singleCell') && params.singleCell
    def singleCellKitArg = params.containsKey('singleCellKit') && params.singleCellKit ? "--single-cell-kit ${params.singleCellKit}" : ''
    def seqspecMd5Enabled = !params.containsKey('seqspecMd5') || params.seqspecMd5
    def seqspecTemplateArg = params.containsKey('seqspecTemplate') && params.seqspecTemplate ? "--template ${params.seqspecTemplate}" : ''
    def seqspecVariablesArg = params.containsKey('seqspecVariables') && params.seqspecVariables ? "--variables ${params.seqspecVariables}" : ''
    """
    . ${params.scriptEnv}
    python ${scripts}/generate_seqspec.py \
        --template-dir ${templates}/seqspec \
        --read-type ${params.readType} \
        ${singleCellEnabled ? '--single-cell' : ''} \
        ${singleCellKitArg} \
        ${seqspecMd5Enabled ? '' : '--no-md5'} \
        ${seqspecTemplateArg} \
        ${seqspecVariablesArg} \
        --fastq ${inputFastq} \
        --output ${outputSpec}
    """
}

process splitcodeTask {
    tag "${params.sample} splitcode"
    container 'ghcr.io/mortazavilab/dogme-pipeline:latest'
    input:
    tuple path(seqspecFile), path(inputFastq)
    output:
    tuple path("${params.sample}_cDNA.fastq.gz"), path("${params.sample}_umi.fastq.gz"), path("${params.sample}_barcode.fastq.gz"), emit: fastqs
    path "${params.sample}_splitcode_qc.tsv", emit: qc
    path "${params.sample}_splitcode.log", emit: splitcode_log
    publishDir "${params.fastqDir}/single-cell", mode: 'copy'
    script:
    def outputFastq = "${params.sample}.splitcode.fastq.gz"
    """
    . ${params.scriptEnv}
    set -o pipefail
    cp ${projectDir}/templates/splitcode/r1_R.txt .
    cp ${projectDir}/templates/splitcode/r1_T.txt .
    cp ${projectDir}/templates/splitcode/r2_3.txt .
    seqspec index -m rna -s file -t splitcode ${seqspecFile} > ONT.config
    sed -i 's/3:3:3/1:1:1/g' ONT.config
    splitcode -c ONT.config -t 2 ${inputFastq} -o ${outputFastq} 2>&1 | tee "${params.sample}_splitcode.log"
    python ${projectDir}/scripts/process_splitcode_fastqs.py \
        --sample "${params.sample}" \
        --qc-output "${params.sample}_splitcode_qc.tsv"
    gunzip -c "${params.sample}_cDNA.fastq.gz" > _cDNA.fastq
    gunzip -c "${params.sample}_umi.fastq.gz" > _umi.fastq
    gunzip -c "${params.sample}_barcode.fastq.gz" > _barcode.fastq
    splitcode -c ${projectDir}/templates/splitcode/config-correct.txt \
        --nFastqs=2 --gzip \
        -o "${params.sample}_cDNA.fastq.gz,_barcode.from-cdna.fastq.gz" \
        _cDNA.fastq _barcode.fastq -t 2
    mv barcode.fastq.gz "${params.sample}_barcode.extracted.fastq.gz"
    splitcode -c ${projectDir}/templates/splitcode/config-correct.txt \
        --nFastqs=2 --gzip \
        -o "${params.sample}_umi.fastq.gz,_barcode.from-umi.fastq.gz" \
        _umi.fastq _barcode.fastq -t 2
    rm -f barcode.fastq.gz
    splitcode -c ${projectDir}/templates/splitcode/config.mergeRT \
        -o "${params.sample}_barcode.fastq" "${params.sample}_barcode.extracted.fastq.gz" -t 2
    gzip -f "${params.sample}_barcode.fastq"
    rm -f "${params.sample}_barcode.extracted.fastq.gz" _barcode.from-cdna.fastq.gz \
        _barcode.from-umi.fastq.gz
    """
}

process singleCellKallistoTask {
    tag "${genomeName} single-cell"
    input:
    tuple path(cDNAFile), path(umiFile), path(barcodeFile), path(indexFile), path(t2gFile), val(genomeName)
    output:
    path "${params.sample}_${genomeName}"
    publishDir "${params.kallistoDir}/${genomeName}/single-cell", mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    output_dir="${params.sample}_${genomeName}"
    mkdir -p "\${output_dir}"
    kallisto bus --long --threshold 0.8 -x '2,0,24:1,0,10:0,0,0' \\
        -i ${indexFile} -t ${task.cpus} -o "\${output_dir}" \\
        "${cDNAFile}" "${umiFile}" "${barcodeFile}"
    bustools whitelist -o "\${output_dir}/whitelist.txt" "\${output_dir}/output.bus"
    bustools correct -w "\${output_dir}/whitelist.txt" \\
        -o "\${output_dir}/corrected.bus" "\${output_dir}/output.bus"
    bustools sort -t ${task.cpus} "\${output_dir}/corrected.bus" \\
        -o "\${output_dir}/sorted.bus"
    bustools count "\${output_dir}/sorted.bus" \\
        -t "\${output_dir}/transcripts.txt" \\
        -e "\${output_dir}/matrix.ec" \\
        -o "\${output_dir}/count" --cm -m -g ${t2gFile}
    """
}

process makeKallistoRefsTask {
    tag "${genomeName}"
    input:
    tuple val(genomeName), val(genomeFasta), val(genomeGtf)
    output:
    tuple val(genomeName), path("${genomeName}.cdna.fa"), path("${genomeName}.introns.fa"), path("${genomeName}.t2g"), val(genomeFasta)
    publishDir "${params.kallistoDir}/refs", mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    python ${projectDir}/scripts/makeKallistoRefs.py --name ${genomeName} --fasta ${genomeFasta} --gtf ${genomeGtf}
    """
}

process kallistoIndexTask {
    tag "${genomeName}"
    input:
    tuple val(genomeName), path(cdnaFa), path(intronsFa), path(t2gFile), val(genomeFasta)
    output:
    tuple val(genomeName), path("${genomeName}.idx"), path(t2gFile)
    publishDir "${params.kallistoDir}/refs", mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    kallisto index -t 8 -i ${genomeName}.idx -k 63 ${cdnaFa}
    """
}

process kallistoTask {
    tag "${genomeName}"
    input:
    tuple path(inputFile), path(indexFile), path(t2gFile), val(genomeName)
    output:
    path "${params.sample}_${genomeName}"
    publishDir "${params.kallistoDir}/${genomeName}", mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    mkdir -p ${params.sample}_${genomeName}
    kallisto bus --long --threshold 0.8 -x bulk -i ${indexFile} -t ${task.cpus} -o ${params.sample}_${genomeName} "${inputFile}"
    bustools sort -t ${task.cpus} ${params.sample}_${genomeName}/output.bus -o ${params.sample}_${genomeName}/sorted.bus
    bustools count ${params.sample}_${genomeName}/sorted.bus -t ${params.sample}_${genomeName}/transcripts.txt  -e ${params.sample}_${genomeName}/matrix.ec  -o ${params.sample}_${genomeName}/count --cm -m -g ${t2gFile}
    kallisto quant-tcc -t ${task.cpus} --long -P ONT ${params.sample}_${genomeName}/count.mtx -i ${indexFile} -f ${params.sample}_${genomeName}/flens.txt -e ${params.sample}_${genomeName}/count.ec.txt -o ${params.sample}_${genomeName}
    """
}

process demuxKallistoTask {
    tag "${sampleName}.${genomeName}"
    input:
    tuple path(inputFile), path(indexFile), path(t2gFile), val(genomeName), val(sampleName)
    output:
    path "${sampleName}_${genomeName}"
    publishDir "${params.kallistoDir}/${genomeName}", mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    mkdir -p ${sampleName}_${genomeName}
    kallisto bus --long --threshold 0.8 -x bulk -i ${indexFile} -t ${task.cpus} -o ${sampleName}_${genomeName} "${inputFile}"
    bustools sort -t ${task.cpus} ${sampleName}_${genomeName}/output.bus -o ${sampleName}_${genomeName}/sorted.bus
    bustools count ${sampleName}_${genomeName}/sorted.bus -t ${sampleName}_${genomeName}/transcripts.txt -e ${sampleName}_${genomeName}/matrix.ec -o ${sampleName}_${genomeName}/count --cm -m -g ${t2gFile}
    kallisto quant-tcc -t ${task.cpus} --long -P ONT ${sampleName}_${genomeName}/count.mtx -i ${indexFile} -f ${sampleName}_${genomeName}/flens.txt -e ${sampleName}_${genomeName}/count.ec.txt -o ${sampleName}_${genomeName}
    """
}

// The splitmodification task generates bed files for each DNA modification. 
// Modifications are identified by letters: 5mCG (m), 6mA (a),and hydroxymethylation (h). 
// The files generated by modkit are grepped for the letter codes. 
process splitModificationTask {
    input:
    path inputFile
    output:
    path "*.filtered*.bed.gz", optional: true
    publishDir params.bedDir, mode: 'copy'
    script:
    """
    . ${params.scriptEnv}
    input_prefix="\$(basename "${inputFile}" .bed.gz)"
    if [[ "${params.readType}" == "DNA" ]]; then 
        # Extract 5mCG (methylation)
        if gzip -cd "${inputFile}" | grep -q -w 'm'; then
            gzip -cd "${inputFile}" | grep -w 'm' | gzip -c > "\${input_prefix}.5mCG.filtered.bed.gz"
        fi
        # Extract 5hmCG (hydroxymethylation)
        if gzip -cd "${inputFile}" | grep -q -w 'h'; then
            gzip -cd "${inputFile}" | grep -w 'h' | gzip -c > "\${input_prefix}.5hmCG.filtered.bed.gz"
        fi
        # Extract 6mA
        if gzip -cd "${inputFile}" | grep -q -w 'a'; then
            gzip -cd "${inputFile}" | grep -w 'a' | gzip -c > "\${input_prefix}.6mA.filtered.bed.gz"
        fi
    elif [[ "${params.readType}" == "RNA" ]]; then
        base_name="\${input_prefix}"
        # Extract m6A modifications (Plus & Minus strands)
        if gzip -cd "${inputFile}" | grep -q -w 'a'; then
            gzip -cd "${inputFile}" | grep -w 'a' | gzip -c > "\${base_name/filtered*/m6A.filtered}.bed.gz"
        fi
        # Extract inosine modifications (Plus & Minus strands)
        if gzip -cd "${inputFile}" | grep -q -w '17596'; then
            gzip -cd "${inputFile}" | grep -w '17596' | gzip -c > "\${base_name/filtered*/inosine.filtered}.bed.gz"
        fi
        # Extract pseudouridine (pseU) modifications (Plus & Minus strands)
        if gzip -cd "${inputFile}" | grep -q -w '17802'; then
            gzip -cd "${inputFile}" | grep -w '17802' | gzip -c > "\${base_name/filtered*/pseU.filtered}.bed.gz"
        fi
        # Extract m5C modifications (Plus & Minus strands)
        if gzip -cd "${inputFile}" | grep -q -w 'm'; then
            gzip -cd "${inputFile}" | grep -w 'm' | gzip -c > "\${base_name/filtered*/m5C.filtered}.bed.gz"
        fi
        # Extract Nm modifications (Plus & Minus strands)
        if gzip -cd "${inputFile}" | grep -qE -w '19228|19229|19227|69426'; then
            gzip -cd "${inputFile}" | grep -Ew '19228|19229|19227|69426' | gzip -c > "\${base_name/filtered*/Nm.filtered}.bed.gz"
        fi
    fi
    """
}

process generateReport {
    tag "Generate metadata and QC reports for ${params.sample}"

    publishDir params.topDir, mode: 'copy'

input:
    val reportInputDir
    val completion
    output:
    path "inventory_report.tsv", emit: inventory_report
    path "qc_summary.csv",       emit: qc_report

    script:
    """
    python ${projectDir}/scripts/generate_report.py \\
        --input_dir ${reportInputDir} \
        --output_inventory inventory_report.tsv \\
        --output_qc qc_summary.csv \\
        --sample ${params.sample}
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

process gtfToJunctionBed {
    input:
    path gtf_file
    output:
    path "*.junctions.bed"
    script:
    """
    python ${projectDir}/scripts/gtf_to_junction_bed.py ${gtf_file} > ${gtf_file.simpleName}.junctions.bed
    """
}

process annotateRNATask {
    input:
    tuple path(bam), path(bai), val(genomeName), path(gtf)
    output:
    val true, emit: completion
    path "*.annotated.ba*"
    path "*.csv"
    path "*_dogme*"
    path "*.log"
    publishDir params.annotDir, mode: 'copy'
    script:
    def outputPrefix = bam.baseName
    """
    . ${params.scriptEnv}
    # If pipeline is running with CDNA read type, pass -CDNA to annotateRNA
    if [[ "${params.readType}" == "CDNA" ]]; then
        cdna_opt="-CDNA"
    else
        cdna_opt=""
    fi

    python ${projectDir}/scripts/annotateRNA.py \
        --bam ${bam} \
        --gtf ${gtf} \
        --out ${outputPrefix} \
        --threads ${task.cpus} \$cdna_opt \
        --novel_prefix "${outputPrefix}" 2> ${outputPrefix}.log
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

    if (model_name.contains('6mA')) {
        reportCompletion = splitResults.collect().combine(consolidatedBeds.collect())
    } else {
        reportCompletion = splitResults.collect()
    }

    emit:
    completion = reportCompletion
}

workflow demuxModificationWorkflow {
    take:
    mapped_bams_ch

    main:
    if (params.readType == 'RNA') {
        mappedBamsForStrands = mapped_bams_ch.map { bam, bai, genomeName, sampleName ->
            tuple(bam, genomeName, sampleName)
        }
        strands = demuxSeparateStrandsTask(mappedBamsForStrands)
        combinedStrand = strands.plus_strand.concat(strands.minus_strand)
            .map { bam, bai, genomeName, sampleName -> tuple(bam, bai, genomeName, sampleName) }
        bedfiles = demuxModkitTask(combinedStrand)
    } else if (params.readType == 'DNA') {
        bedfiles = demuxModkitTask(mapped_bams_ch)
    }

    filterbeds = filterbedTask(bedfiles)
    splitResults = splitModificationTask(filterbeds)

    emit:
    completion = splitResults.collect()
}

workflow kallistoWorkflow {
    take:
    unmapped_bams_ch

    main:
    def singleCellEnabled = params.containsKey('singleCell') && params.singleCell
    fastqFile = extractfastqTask(unmapped_bams_ch)

    seqspecFile = generateSeqspecTask(
        fastqFile,
        file("${projectDir}/templates"),
        file("${projectDir}/scripts"),
    )

    if (params.readType == 'CDNA' && singleCellEnabled) {
        splitcodeInput = fastqFile.combine(seqspecFile)
            .map { fastq, spec -> tuple(spec, fastq) }
        splitcodeTask(splitcodeInput)
    }

    if (!params.kallistoIndex || !params.t2g) {
        def kallisto_refs_ch = nextflow.Channel.fromList(params.genome_annot_refs)
            .map { ref -> tuple(ref.name, ref.genome, ref.annot) }
        refFiles = makeKallistoRefsTask(kallisto_refs_ch)
        indexFiles = kallistoIndexTask(refFiles)
        if (singleCellEnabled) {
            singleCellInput = splitcodeTask.out.fastqs.combine(indexFiles)
                .map { cDNA, umi, barcode, genomeName, idx, t2g ->
                    tuple(cDNA, umi, barcode, idx, t2g, genomeName)
                }
            terminalKallisto = singleCellKallistoTask(singleCellInput)
        } else {
            kallistoInput = fastqFile.combine(indexFiles)
                .map { fastq, genomeName, idx, t2g -> tuple(fastq, idx, t2g, genomeName) }
            terminalKallisto = kallistoTask(kallistoInput)
        }
    } else {
        if (singleCellEnabled) {
            singleCellInput = splitcodeTask.out.fastqs.map { cDNA, umi, barcode ->
                tuple(cDNA, umi, barcode, file(params.kallistoIndex), file(params.t2g), 'prebuilt')
            }
            terminalKallisto = singleCellKallistoTask(singleCellInput)
        } else {
            kallistoInput = fastqFile.map { fastq ->
                tuple(fastq, file(params.kallistoIndex), file(params.t2g), 'prebuilt')
            }
            terminalKallisto = kallistoTask(kallistoInput)
        }
    }

    emit:
    completion = terminalKallisto.collect().combine(seqspecFile.collect())
}

workflow kallistoFastqWorkflow {
    take:
    fastq_ch

    main:
    fastqFile = fastq_ch

    if (!params.kallistoIndex || !params.t2g) {
        def kallisto_refs_ch = nextflow.Channel.fromList(params.genome_annot_refs)
            .map { ref -> tuple(ref.name, ref.genome, ref.annot) }
        refFiles = makeKallistoRefsTask(kallisto_refs_ch)
        indexFiles = kallistoIndexTask(refFiles)
        kallistoInput = fastqFile.combine(indexFiles)
            .map { fastq, genomeName, idx, t2g -> tuple(fastq, idx, t2g, genomeName) }
        terminalKallisto = kallistoTask(kallistoInput)
    } else {
        kallistoInput = fastqFile.map { fastq ->
            tuple(fastq, file(params.kallistoIndex), file(params.t2g), 'prebuilt')
        }
        terminalKallisto = kallistoTask(kallistoInput)
    }

    emit:
    completion = terminalKallisto.collect()
}

workflow demuxKallistoWorkflow {
    take:
    fastq_ch

    main:
    if (!params.kallistoIndex || !params.t2g) {
        def kallisto_refs_ch = nextflow.Channel.fromList(params.genome_annot_refs)
            .map { ref -> tuple(ref.name, ref.genome, ref.annot) }
        refFiles = makeKallistoRefsTask(kallisto_refs_ch)
        indexFiles = kallistoIndexTask(refFiles)
        kallistoInput = fastq_ch.combine(indexFiles)
            .map { fastq, sampleName, genomeName, idx, t2g ->
                tuple(fastq, idx, t2g, genomeName, sampleName)
            }
        terminalKallisto = demuxKallistoTask(kallistoInput)
    } else {
        kallistoInput = fastq_ch.map { fastq, sampleName ->
            tuple(fastq, file(params.kallistoIndex), file(params.t2g), 'prebuilt', sampleName)
        }
        terminalKallisto = demuxKallistoTask(kallistoInput)
    }

    emit:
    completion = terminalKallisto.collect()
}

workflow mainWorkflow {
    take:
    theVersion
    theModel 
    modelDirectory
    
    main: 
    // Only run dorado download if the model directory does not already exist
    def modelPath
    if (!new File(modelDirectory).exists()) {
        modelPath = doradoDownloadTask(modelDirectory, theModel)
    } else {
        modelPath = nextflow.Channel.value(modelDirectory)
    }
    softwareVTask(theVersion, modelPath)
    def pod5_files_ch = nextflow.Channel.fromPath("${params.podDir}/*.pod5")
    bamFiles = doradoTask(pod5_files_ch, modelDirectory, modelPath, theModel).collectFile()
    fileCount = bamFiles.map { it.size() }.first()
    unmappedbam = mergeBamsTask(fileCount)
    def genome_annot_ch = nextflow.Channel.fromList(params.genome_annot_refs)
    def demuxAnnotationBams = Channel.empty()
    def analysisCompletions = Channel.empty()
    if (params.kitName) {
        demuxBamFiles = doradoDemuxTask(unmappedbam).demuxed_bams
            .flatten()
        demuxSummary = demuxSummaryTask(demuxBamFiles.collect())
        analysisCompletions = analysisCompletions.mix(demuxSummary.collect())
        classifiedBams = demuxBamFiles
            .filter { bam ->
                def name = bam.simpleName.toLowerCase()
                !name.contains('unclassified') && !name.contains('no_barcode')
            }
        if (params.keepBarcodes != null) {
            classifiedBams = selectTopDemuxBamsTask(classifiedBams.collect()).selected_bams
                .flatten()
                .map { bam -> tuple(bam, bam.baseName) }
        } else {
            classifiedBams = demuxBamFiles.map { bam ->
                def barcodeDirectory = bam.parent.name
                def barcodeMatch = barcodeDirectory =~ /(?i)^barcode(\d+)$/
                if (!barcodeMatch.matches()) {
                    throw new IllegalArgumentException("Expected a Dorado barcode directory for '${bam}', found '${barcodeDirectory}'")
                }
                tuple(bam, "bc${barcodeMatch[0][1]}")
            }
        }
        prepareDemuxBamTask(classifiedBams)
        classifiedUnmapped = prepareDemuxBamTask.out.classified_unmapped
            .map { bam ->
                def sampleName = bam.name.replaceFirst(/\.unmapped\.bam$/, '')
                tuple(bam, sampleName)
            }
        if (params.readType == 'RNA' || params.readType == 'CDNA') {
            demuxFastq = demuxExtractfastqTask(classifiedUnmapped)
            demuxSeqspec = demuxGenerateSeqspecTask(
                demuxFastq,
                file("${projectDir}/templates"),
                file("${projectDir}/scripts"),
            )
            demuxKallisto = demuxKallistoWorkflow(demuxFastq)
            analysisCompletions = analysisCompletions.mix(demuxSeqspec.collect())
            analysisCompletions = analysisCompletions.mix(demuxKallisto.completion)
        }
        demuxMappedBams = classifiedUnmapped
            .combine(genome_annot_ch)
            .map { bam, sampleName, ref ->
                tuple(bam, sampleName, ref.genome, ref.annot, ref.name)
            }
        demuxMappedResults = demuxMinimapTask(demuxMappedBams)
        if (params.readType == 'RNA' || params.readType == 'CDNA') {
            demuxAnnotationBams = demuxMappedResults
                .map { bam, bai, genomeName, sampleName -> tuple(bam, bai, genomeName) }
        }
        if (params.readType == 'RNA' || params.readType == 'DNA') {
            demuxModifications = demuxModificationWorkflow(demuxMappedResults)
            analysisCompletions = analysisCompletions.mix(demuxModifications.completion)
        }
    }
    
    if (params.readType == 'RNA' || params.readType == 'CDNA') {
        kallistoResults = kallistoWorkflow(unmappedbam)
        analysisCompletions = analysisCompletions.mix(kallistoResults.completion)
    }

    unmappedBams = unmappedbam.combine(genome_annot_ch).map { bam, ref ->
        tuple(bam, ref.genome, ref.annot, ref.name)
    }
    minimapTask(unmappedBams)
    def mappedBams = minimapTask.out.mapped_bams

    if (params.readType == 'RNA' || params.readType == 'DNA') {
        modifications = modificationWorkflow(mappedBams, theModel)
        analysisCompletions = analysisCompletions.mix(modifications.completion)
    } else {
        analysisCompletions = analysisCompletions.mix(mappedBams.collect())
    }

    if (params.readType == 'RNA' || params.readType == 'CDNA') {
        annotations = annotateRNAWorkflow(mappedBams.mix(demuxAnnotationBams))
        analysisCompletions = analysisCompletions.mix(annotations.completion)
    }

    generateReport(params.topDir, analysisCompletions.collect())
}

workflow basecallWorkflow {
    take:
    theVersion
    theModel 
    modelDirectory
    
    main: 
    modelPath = doradoDownloadTask(modelDirectory, theModel)
    softwareVTask(theVersion, modelPath)
    def pod5_files_ch = nextflow.Channel.fromPath("${params.podDir}/*.pod5")
    bamFiles = doradoTask(pod5_files_ch, modelDirectory, modelPath, theModel).collectFile()
    fileCount = bamFiles.map { it.size() }.first()
    unmappedbam = mergeBamsTask(fileCount)
    if (params.kitName) {
        demuxedBams = doradoDemuxTask(unmappedbam)
        classifiedBams = demuxedBams.demuxed_bams
            .flatten()
            .filter { bam ->
                def name = bam.simpleName.toLowerCase()
                !name.contains('unclassified') && !name.contains('no_barcode')
            }
            .map { bam -> tuple(bam, bam.simpleName) }
        prepareDemuxBamTask(classifiedBams)
    }
}

workflow remapWorkflow {
    take:
    theVersion
    theModel 
    modelDirectory

    main:
    def unmappedbam = nextflow.Channel.fromPath("${params.bamDir}/*.unmapped.bam")
    def genomeAnnotChannel = nextflow.Channel.fromList(params.genome_annot_refs)
    unmappedBams = unmappedbam.combine(genomeAnnotChannel).map { bam, ref ->
        tuple(bam, ref.genome, ref.annot, ref.name)
    }
    minimapTask(unmappedBams)
    def mappedBams = minimapTask.out.mapped_bams
    def analysisCompletions = Channel.empty()

    if (params.readType == 'RNA' || params.readType == 'DNA') {
        modifications = modificationWorkflow(mappedBams, theModel)
        analysisCompletions = analysisCompletions.mix(modifications.completion)
    } else {
        analysisCompletions = analysisCompletions.mix(mappedBams.collect())
    }

    if (params.readType == 'RNA' || params.readType == 'CDNA') {
        kallistoResults = kallistoWorkflow(unmappedbam)
        annotations = annotateRNAWorkflow(mappedBams)
        analysisCompletions = analysisCompletions.mix(kallistoResults.completion)
        analysisCompletions = analysisCompletions.mix(annotations.completion)
    }

    generateReport(params.topDir, analysisCompletions.collect())
}

workflow fastqCDNAWorkflow {
    take:
    theVersion
    theModel
    modelDirectory

    main:
    def fastqCandidates = [
        new File(params.fastqDir, "${params.sample}.fastq.gz"),
        new File(params.fastqDir, "${params.sample}.fastq")
    ].findAll { it.exists() }

    if (!fastqCandidates) {
        throw new IllegalArgumentException("No FASTQ found for sample '${params.sample}' in ${params.fastqDir}. Expected ${params.sample}.fastq.gz or ${params.sample}.fastq")
    }

    if (fastqCandidates.size() > 1) {
        def matches = fastqCandidates.collect { it.name }.join(', ')
        throw new IllegalArgumentException("Multiple FASTQ inputs found for sample '${params.sample}' in ${params.fastqDir}: ${matches}")
    }

    def inputFastq = nextflow.Channel.fromPath(fastqCandidates[0].path, checkIfExists: true)
    def unmappedbam = samtoolsImportTask(inputFastq)

    kallistoResults = kallistoFastqWorkflow(inputFastq)

    def genome_annot_ch = nextflow.Channel.fromList(params.genome_annot_refs)

    unmappedBams = unmappedbam.combine(genome_annot_ch).map { bam, ref ->
        tuple(bam, ref.genome, ref.annot, ref.name)
    }
    minimapTask(unmappedBams)
    def mappedBams = minimapTask.out.mapped_bams

    annotations = annotateRNAWorkflow(mappedBams)
    generateReport(params.topDir, kallistoResults.completion.mix(annotations.completion).collect())
}

workflow reportsWorkflow {
    take:
    theVersion
    modelDirectory
    
    main:
    softwareVTask(theVersion, modelDirectory)
    
    generateReport(params.topDir, nextflow.Channel.of(true))
}

workflow annotateRNAWorkflow {
    take:
    mapped_bams_ch

    main:
    // 1. Create the GTF channel with a new name: 'annotation_gtf_ch'
    def annotation_gtf_ch = Channel
        .fromList(params.genome_annot_refs)
        .map { ref -> tuple(ref.name, file(ref.annot)) }

    // 2. Prepare the BAM channel for grouping
    def bams_for_grouping = mapped_bams_ch
        .map { bam, bai, genomeName -> tuple(genomeName, bam, bai) }

    // 3. Group all BAMs by their genome name
    def grouped_bams_ch = bams_for_grouping.groupTuple()

    // 4. Combine with the renamed channel: 'annotation_gtf_ch'
    def combined_ch = grouped_bams_ch.combine(annotation_gtf_ch, by: 0)

    // 5. "Un-roll" the grouped structure back into a flat channel of 4 items
    def mappedBamsWithGtf = combined_ch.flatMap { genomeName, bams, bais, gtf_file ->
        def results = []
        // Loop through each BAM in the group
        for( int i = 0; i < bams.size(); i++ ) {
            // Create a final tuple for each BAM: (bam, bai, genomeName, gtf_file)
            results.add( tuple(bams[i], bais[i], genomeName, gtf_file) )
        }
        return results
    }

    annotationResults = annotateRNATask(mappedBamsWithGtf)

    emit:
    completion = annotationResults.completion.collect()
}