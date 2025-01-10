#!/usr/bin/env nextflow


nextflow.enable.dsl=2

include { modWorkflow } from './nanoporeModule'


workflow {
    modWorkflow('sup,inosine_m6A,pseU,m5C')
}



// Define other processes similarly, ensuring they take `rnaMod` as an input
workflow.onComplete {
    // Execute cleanup command
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

