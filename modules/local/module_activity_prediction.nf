process MODULE_ACTIVITY_PREDICTION {
    tag  'module_activity_prediction'
    tag  "${if (workflow.stubRun) 'module_activity_prediction_stub' else 'module_activity_prediction'}"
    label 'process_high'
    container 'registry.gitlab.com/bzhanglab/python:3.8.13'

    input:
    path config_file
    path data_tgz_file
    path module_info_file
    path mutation_file
    path network_edge_list

    output:
    path '*.pdf', emit: module_activity_plots
    path '*.{tsv,json}', emit: module_activity_data
    path 'selected_modules.txt', emit: selected_modules
    path 'versions.yml', emit: versions

    when:
    params.mutation_file && params.module_activity_dataset_name

    script:
    """
    extract_data_file.py \\
        --dataset-name ${params.module_activity_dataset_name} \\
        --config-file ${config_file} \\
        --data-file ${data_tgz_file} \\
        --output-file data.tsv
    module_activity_prediction.py \\
        --dataset-name ${params.module_activity_dataset_name} \\
        --module-info ${module_info_file} \\
        --data-file data.tsv \\
        --mutation-file ${mutation_file}
    rm -rf data.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: 3.8.13
    END_VERSIONS
    """

    stub:
    """
    wget "https://drive.google.com/uc?id=1qf8YvUM-d12eU2zg7Ek5Stg41HwsN2K-&export=download&confirm=9iBg" -O module_activity_prediction_results.tgz
    tar -xzf module_activity_prediction_results.tgz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: 3.8.13
    END_VERSIONS
    """
}
