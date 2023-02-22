#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

WorkflowMain.initialise(workflow, params, log)

include { IFUNMAP } from './workflows/ifunmap'

workflow MOCORE_IFUNMAP {
    IFUNMAP ()
}

workflow {
    MOCORE_IFUNMAP ()
}
