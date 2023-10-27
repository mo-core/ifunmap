#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

WorkflowMain.initialise(workflow, params, log)

include { CHECK_CYTOSCAPE} from './modules/local/check_cytoscape'
include { IFUNMAP } from './workflows/ifunmap'

workflow MOCORE_IFUNMAP {
    IFUNMAP ()
}

workflow {
    exitCode = CHECK_CYTOSCAPE()
    if (exitCode == 0) {
        MOCORE_IFUNMAP ()
    } else {
        // print error message
        println "Cytoscape is not running..."
    }
}
