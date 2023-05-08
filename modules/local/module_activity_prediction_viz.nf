process MODULE_ACTIVITY_PREDICTION_VIZ {
    tag  'module_activity_prediction_viz'
    tag  "${if (workflow.stubRun) 'module_activity_prediction_viz_stub' else 'module_activity_prediction_viz'}"
    executor 'local'
    container 'registry.gitlab.com/bzhanglab/python:3.8.13'
    containerOptions '--net host --privileged'

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

   stub:
    """
    wget "https://drive.google.com/uc?id=1WPJpQguwEUK6dYPGICWPZ5q-HV1UAJBO&export=download&confirm=9iBg" -O module_activity_prediction_viz_results.tgz
    tar -xzf module_activity_prediction_viz_results.tgz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: 3.8.13
    END_VERSIONS
    """

}
