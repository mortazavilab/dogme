#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process openChromatinTask {
    tag "${genomeName} ${chromosome} ${outputType}"

    input:
    tuple val(chromosome), val(size), path(bam_file), path(bai_file), val(genomeName), val(outputType)

    output:
    path "${genomeName}/${chromosome}/${chromosome}*"

    publishDir "${params.topDir}/open_chromatin", mode: 'copy'

    script:
    """
    #!/bin/bash
    . ${params.scriptEnv}

    target_chromosome="${chromosome}"
    output_dir="${genomeName}/\${target_chromosome}"
    mkdir -p "\$output_dir"

    chunk_size=1000000

    if [ -z "${size}" ] || [ "${size}" -eq 0 ]; then
        echo "Warning: Chromosome '\$target_chromosome' has a size of 0. Skipping."
        exit 0
    fi

    final_file="\${output_dir}/\${target_chromosome}.bg"
    final_region_file="\${output_dir}/\${target_chromosome}.bed"

    chunk_count=0

    if [[ "${outputType}" == "bg" ]]; then
        echo -e "#chrom\\tstart\\tstop\\tprob" > "\$final_file"
    elif [[ "${outputType}" == "bed" ]]; then
        > "\$final_region_file"
    fi

    for ((start=1; start<=${size}; start+=chunk_size)); do
        end=\$((start + chunk_size - 1))
        if [ "\$end" -gt "${size}" ]; then
            end="${size}"
        fi

        ((chunk_count++))
        chunk="\${target_chromosome}:\${start}-\${end}"
        suffix=\$(printf "%04d" "\$chunk_count")

        if [[ "${outputType}" == "bg" ]]; then
            log_file="\${output_dir}/temp_\${target_chromosome}_\${start}_\${end}_\${suffix}.log"
            output_file="\${output_dir}/temp_\${target_chromosome}_\${start}_\${end}_\${suffix}.bg"
            echo "Predicting signal for \$chunk, output: \$output_file"
            modkit open-chromatin predict "${bam_file}" --model "\${MODKITMODEL}" --log "\$log_file" --no-header -o "\$output_file" --region "\$chunk"
            cat "\$output_file" >> "\$final_file"
        elif [[ "${outputType}" == "bed" ]]; then
            region_file="\${output_dir}/temp_\${target_chromosome}_\${start}_\${end}_\${suffix}.bed"
            echo "Predicting elements for \$chunk, output: \$region_file"
            modkit open-chromatin predict "${bam_file}" --model "\${MODKITMODEL}" --threshold 0.8 -o stdout --region "\$chunk" | bedtools merge -i - > "\$region_file"
            cat "\$region_file" >> "\$final_region_file"
        fi

        if [ "\$end" -eq "${size}" ]; then
            break
        fi
    done
    """
}
