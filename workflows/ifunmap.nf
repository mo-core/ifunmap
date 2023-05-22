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
include { FUNMAP_QC } from '../modules/local/funmap_qc'
include { FUNMAP_RUN } from '../modules/local/funmap_run'
include { ICE } from '../modules/local/ice'
include { CLIQUE_ENRICH } from '../modules/local/clique_enrich'
include { CLIQUE_VIZ } from '../modules/local/clique_viz'
include { NETWORK_ANALYSIS } from '../modules/local/network_analysis'
include { MODULE_ACTIVITY} from '../modules/local/module_activity'
include { MODULE_ACTIVITY_PREDICTION} from '../modules/local/module_activity_prediction'
include { MODULE_ACTIVITY_PREDICTION_VIZ} from '../modules/local/module_activity_prediction_viz'
include { PLOT_DARK_GENES_HEATMAP } from '../modules/local/plot_dark_genes_heatmap'
include { PLOT_DARK_GENE_ENRICH_PIE } from '../modules/local/plot_dark_gene_enrich_pie'
include { DARK_GENE_ENRICH } from '../modules/local/dark_gene_enrich'

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

    if (params.qc_only) {
        FUNMAP_QC(
            INPUT_CHECK.out.config_file,
            ch_data)
        ch_versions = ch_versions.mix(FUNMAP_QC.out.versions)
    } else {
        FUNMAP_RUN (
            INPUT_CHECK.out.config_file,
            ch_data
        )
        ch_versions = ch_versions.mix(FUNMAP_RUN.out.versions)

        // dark gene analysis is optional, default is true
        if (params.dark_gene_tgi) {
            ch_dark_genes = file(params.dark_gene_tgi)
        } else {
            ch_dark_genes = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)
        }
        ch_pubmed_count = file(params.gene_pubmed_count)

        PLOT_DARK_GENES_HEATMAP(
            FUNMAP_RUN.out.funmap_el,
            ch_dark_genes,
            ch_pubmed_count
        )
        ch_versions = ch_versions.mix(PLOT_DARK_GENES_HEATMAP.out.versions)

        DARK_GENE_ENRICH(
            FUNMAP_RUN.out.funmap_el,
            ch_pubmed_count
        )
        ch_versions = ch_versions.mix(DARK_GENE_ENRICH.out.versions)

        PLOT_DARK_GENE_ENRICH_PIE(
            DARK_GENE_ENRICH.out.enrich_results,
            DARK_GENE_ENRICH.out.top_neighbors
        )
        ch_versions = ch_versions.mix(PLOT_DARK_GENE_ENRICH_PIE.out.versions)

        ICE (
            FUNMAP_RUN.out.funmap_el
        )
        ch_versions = ch_versions.mix(ICE.out.versions)

        CLIQUE_ENRICH (
            FUNMAP_RUN.out.funmap_el,
            ICE.out.ice_results
        )
        ch_versions = ch_versions.mix(CLIQUE_ENRICH.out.versions)

        NETWORK_ANALYSIS (
            FUNMAP_RUN.out.funmap_el
        )
        ch_versions = ch_versions.mix(NETWORK_ANALYSIS.out.versions)

        CLIQUE_VIZ (
            FUNMAP_RUN.out.funmap_el,
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

        if (params.mutation_file) { ch_mut = file(params.mutation_file) } else { ch_mut = Channel.empty() }
        MODULE_ACTIVITY_PREDICTION (
            INPUT_CHECK.out.config_file,
            ch_data,
            MODULE_ACTIVITY.out.module_info,
            ch_mut,
            FUNMAP_RUN.out.funmap_el
        )
        ch_versions = ch_versions.mix(MODULE_ACTIVITY_PREDICTION.out.versions)

        MODULE_ACTIVITY_PREDICTION_VIZ (
            FUNMAP_RUN.out.funmap_el,
            MODULE_ACTIVITY_PREDICTION.out.selected_modules,
            MODULE_ACTIVITY.out.module_info
        )
        ch_versions = ch_versions.mix(MODULE_ACTIVITY_PREDICTION_VIZ.out.versions)
    }
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
