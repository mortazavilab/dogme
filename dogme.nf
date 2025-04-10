#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { modWorkflow } from './nanoporeModule'

def getParamOrDefault(param, defaultValue) {
    if (param == null || param == 'null' || param == 'undefined' || !param) {
        return defaultValue
    } else {
        return param
    }
}

// Set the default value at the workflow level
def dogmeVersion = "0.90"
def defaultModDir = "${launchDir}/doradoModels"

workflow {
    modDir = getParamOrDefault(params.modDir, defaultModDir)
    params.modDir = modDir

    def modificationsMap = [
        "RNA": 'inosine_m6A,pseU,m5C',
        "DNA": '5mCG_5hmCG,6mA'
    ]

    theModifications = getParamOrDefault(params.modifications, modificationsMap.get(params.readType, ''))
    theModel = params.accuracy + (theModifications ? ",${theModifications}" : "")

    results = modWorkflow(dogmeVersion, theModel, modDir)
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