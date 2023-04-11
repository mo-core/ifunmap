/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowIfunmap.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.config_file, params.data_file ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.config_file) { ch_config = file(params.config_file) } else { exit 1, 'config file not specified!' }
if (params.data_file) { ch_data = file(params.data_file) } else { exit 1, 'config file not specified!' }

include { INPUT_CHECK } from '../subworkflows/local/input_check'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { FUNMAP } from '../modules/local/funmap'
include { ICE } from '../modules/local/ice'
include { CLIQUE_ENRICH } from '../modules/local/clique_enrich'
include { CLIQUE_VIZ } from '../modules/local/clique_viz'
include { NETWORK_ANALYSIS } from '../modules/local/network_analysis'
include { MODULE_ACTIVITY} from '../modules/local/module_activity'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary

workflow IFUNMAP {

    ch_versions = Channel.empty()

    INPUT_CHECK (
        ch_config,
        ch_data
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    FUNMAP (
        INPUT_CHECK.out.config_file,
        ch_data
    )
    ch_versions = ch_versions.mix(FUNMAP.out.versions)

    ICE (
        FUNMAP.out.funmap_el
    )
    ch_versions = ch_versions.mix(ICE.out.versions)

    CLIQUE_ENRICH (
        FUNMAP.out.funmap_el,
        ICE.out.ice_results
    )
    ch_versions = ch_versions.mix(CLIQUE_ENRICH.out.versions)

    NETWORK_ANALYSIS (
        FUNMAP.out.funmap_el
    )
    ch_versions = ch_versions.mix(NETWORK_ANALYSIS.out.versions)

    CLIQUE_VIZ (
        FUNMAP.out.funmap_el,
        ICE.out.ice_results,
        CLIQUE_ENRICH.out.clique_enrich_results,
        INPUT_CHECK.out.config_file,
        ch_data
    )
    ch_versions = ch_versions.mix(CLIQUE_VIZ.out.versions)

    if (params.tsi) { ch_tsi = file(params.tsi) } else { ch_tsi = Channel.empty() }
    MODULE_ACTIVITY (
        ICE.out.ice_results,
        NETWORK_ANALYSIS.out.netsam_nsm,
        INPUT_CHECK.out.config_file,
        ch_data,
        ch_tsi
    )
    ch_versions = ch_versions.mix(MODULE_ACTIVITY.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // workflow_summary    = WorkflowIfunmap.paramsSummaryMultiqc(workflow, summary_params)
    // ch_workflow_summary = Channel.value(workflow_summary)

    // methods_description    = WorkflowIfunmap.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    // ch_methods_description = Channel.value(methods_description)

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}
