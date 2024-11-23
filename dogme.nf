#!/usr/bin/env nextflow


nextflow.enable.dsl=2

include { modWorkflow as m6aWorkflow } from './nanoporeModule'
include { modWorkflow as pseUWorkflow} from './nanoporeModule'
include { modWorkflow as m5cWorkflow} from './nanoporeModule'


workflow {
    m6aWorkflow('m6aInosine', 'sup,inosine_m6A')
    pseUWorkflow('pseU', 'sup,pseU')
    m5cWorkflow('m5c', 'sup,m5C')
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

