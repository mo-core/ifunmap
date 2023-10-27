#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

WorkflowMain.initialise(workflow, params, log)

include { IFUNMAP } from './workflows/ifunmap'

workflow MOCORE_IFUNMAP {
    IFUNMAP ()
}

workflow {
    def dockerImageName = "cytoscape/cytoscape-desktop"
    def cmd = ["docker", "ps", "--format", "{{.Image}}"]
    def processBuilder = new ProcessBuilder(cmd)
    def process = processBuilder.start()
    def reader = new BufferedReader(new InputStreamReader(process.getInputStream()))
    def line
    def isRunning = false

    while ((line = reader.readLine()) != null) {
        if (line.contains(dockerImageName)) {
            isRunning = true
            break
        }
    }
    process.waitFor()

    if (isRunning) {
        MOCORE_IFUNMAP ()
    } else {
        println "Cytoscape is not running..."
        println "Please start a container based on cytoscape/cytoscape-desktop before running this workflow"
        println "e.g., docker run -p 1234:1234 -v /data:/data -it cytoscape/cytoscape-desktop:3.9.1"
        println "replace /data with the path to your data directory"
    }
}
