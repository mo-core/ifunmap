process MODULE_ACTIVITY {
    tag  "${if (workflow.stubRun) 'module_activity_stub' else 'module_activity'}"
    label 'process_high'
    container 'registry.gitlab.com/bzhanglab/python:3.8.13'
    containerOptions '--privileged'

    input:
    path ice_clique
    path netsam_nsm
    path config_file
    path data_file
    path tsi_file


    output:
    path 'module_info.tsv', emit: module_info
    path 'module_activity_*.*', emit: module_activity_scores
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
       python: 3.8.13
    END_VERSIONS
    """

    stub:
    """
    wget "https://drive.google.com/uc?id=19KhjLjYFCdiSmniag4P2Xvtrp96txclY&export=download&confirm=9iBg" -O module_activity_results.tgz
    tar -xzf module_activity_results.tgz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: 3.8.13
    END_VERSIONS
    """

}
