process MODULE_ACTIVITY_PREDICTION_VIZ {
    tag  'module_activity_prediction_viz'
    executor 'local'
    container 'registry.gitlab.com/bzhanglab/python:3.8.13'
    containerOptions '--net host'

    input:
    path network_edge_list
    path module_list
    path module_info

    output:
    path '*.pdf', emit: plots
    path '*.cys', emit: cys
    path 'versions.yml', emit: versions

    when:
    params.mutation_file && params.module_activity_dataset_name

    script:
    """
    module_activity_prediction_viz.py \\
        --network-edge-list ${network_edge_list} \\
        --module-list ${module_list} \\
        --module-info ${module_info}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: 3.8.13
    END_VERSIONS
    """

}
