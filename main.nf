#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

WorkflowMain.initialise(workflow, params, log)

include { CHECK_CYTOSCAPE} from './modules/local/check_cytoscape'
include { IFUNMAP } from './workflows/ifunmap'

workflow MOCORE_IFUNMAP {
    IFUNMAP ()
}

workflow {
    if (CHECK_CYTOSCAPE()) {
        MOCORE_IFUNMAP ()
    }
}
