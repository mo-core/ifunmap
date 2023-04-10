process MODULE_ACTIVITY {
    tag  "${if (workflow.stubRun) 'module_activity_stub' else 'module_activity'}"
    label 'process_high'

    container 'registry.gitlab.com/bzhanglab/python:3.7.6'

    input:
    path ice_clique
    path netsam_nsm
    path config_file
    path data_file
    path tsi_file


    output:
    path 'module_dict.npy', emit: module_dict
    path 'module_activity_scores.*', emit: module_activity_scores
    path 'versions.yml', emit: versions

    when:
    params.tsi && params.module_activity_dataset_name

    script:
    """
    extract_data_file.py \\
        --dataset-name ${params.module_activity_dataset_name} \\
        --config-file ${config_file} \\
        --data-file ${data_file} \\
        --output-file data.tsv
    module_activity.py \\
        --ice-clique ${ice_clique} \\
        --netsam-nsm ${netsam_nsm} \\
        --data-file data.tsv \\
        --tsi-file ${tsi_file}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: 3.7.6
    END_VERSIONS
    """

    // stub:
    // """
    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //   python: 3.7.6
    // END_VERSIONS
    // """
}
