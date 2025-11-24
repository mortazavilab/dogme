#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { mainWorkflow } from './nanoporeModule'
include { reportsWorkflow } from './nanoporeModule'
include { remapWorkflow } from './nanoporeModule'
include { basecallWorkflow } from './nanoporeModule'
include { modificationWorkflow } from './nanoporeModule'
include { annotateRNAWorkflow } from './nanoporeModule'

def getParamOrDefault(param, defaultValue) {
    if (param == null || param == 'null' || param == 'undefined' || !param) {
        return defaultValue
    } else {
        return param
    }
}

// Set the default value at the workflow level
def dogmeVersion = "1.2.3"
def defaultModDir = "${launchDir}/doradoModels"

// Define modificationsMap once here, to be reused across workflows
def modificationsMap = [
    "RNA": '2OmeG,m5C_2OmeC,inosine_m6A_2OmeA,pseU_2OmeU',
    "DNA": '5mCG_5hmCG,6mA'
]

workflow {
    modDir = getParamOrDefault(params.modDir, defaultModDir)
    params.modDir = modDir

    theModifications = getParamOrDefault(params.modifications, modificationsMap.get(params.readType, ''))
    theModel = params.accuracy + (theModifications ? ",${theModifications}" : "")

    results = mainWorkflow(dogmeVersion, theModel, modDir)
}

workflow basecall {
    modDir = getParamOrDefault(params.modDir, defaultModDir)
    params.modDir = modDir

    theModifications = getParamOrDefault(params.modifications, modificationsMap.get(params.readType, ''))
    theModel = params.accuracy + (theModifications ? ",${theModifications}" : "")

    results = basecallWorkflow(dogmeVersion, theModel, modDir)
}

workflow remap {
    modDir = getParamOrDefault(params.modDir, defaultModDir)
    params.modDir = modDir

    theModifications = getParamOrDefault(params.modifications, modificationsMap.get(params.readType, ''))
    theModel = params.accuracy + (theModifications ? ",${theModifications}" : "")

    results = remapWorkflow(dogmeVersion, theModel, modDir)
}

workflow modkit {
    theModifications = getParamOrDefault(params.modifications, modificationsMap.get(params.readType, ''))
    theModel = params.accuracy + (theModifications ? ",${theModifications}" : "")

    // 1. Create a channel from all BAM files in the specified directory.
    // 2. Filter out any containing "unmapped".
    // 3. Filter to keep only those that have a corresponding .bai index file.
    // 4. Map the results to the tuple format expected by analysisWorkflow.
    mappedBams = Channel.fromPath("${params.bamDir}/*.bam")
                         .filter { bam -> !bam.name.contains('unmapped') }
                         .filter { bam -> file(bam.toString() + '.bai').exists() }
                         .map { bam ->
                             def bai = file(bam.toString() + '.bai')
                             // Assumes the filename format is ${params.sample}.${genomeName}.bam
                             // This line removes the sample prefix to extract the genome name.
                             def genomeName = bam.baseName.replaceFirst("^${params.sample}\\.", "")
                             return tuple(bam, bai, genomeName)
                         }

    modificationWorkflow(mappedBams, theModel)
}

workflow annotateRNA {
    // 1. Create a channel from all mapped BAM files in the specified directory.
    // 2. Filter out any containing "unmapped".
    // 3. Filter to keep only those that have a corresponding .bai index file.
    // 4. Map the results to the tuple format expected by annotateRNAWorkflow.
    mappedBams = Channel.fromPath("${params.bamDir}/*.bam")
        .filter { bam -> !bam.name.contains('unmapped') }
        .filter { bam -> file(bam.toString() + '.bai').exists() }
        .map { bam ->
            def bai = file(bam.toString() + '.bai')
            // Assume filename format is: sample_name.genome.bam
            // Extract genomeName as the last dot-separated field of the base name
            def baseName = bam.baseName
            def parts = baseName.tokenize('.')
            def genomeName = parts.size() > 1 ? parts[-1] : baseName
            return tuple(bam, bai, genomeName)
        }


    annotateRNAWorkflow(mappedBams)
}

workflow reports {
    reportsWorkflow(dogmeVersion, params.modDir)
}

workflow.onComplete {
    def cleanCommand = "rm -rf ${workflow.workDir}/*"
    println "Cleaning up work directory..."
    def process = cleanCommand.execute()
    process.waitFor()
    if (process.exitValue() == 0) {
        println "Work directory cleaned successfully."
    } else {
        println "Failed to clean work directory."
        println process.err.text
    }
}
