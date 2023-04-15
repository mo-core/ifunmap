process MODULE_ACTIVITY_PREDICTION {
    tag  'module_activity_prediction'
    label 'process_high'
    container 'registry.gitlab.com/bzhanglab/python:3.7.6'

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
       python: 3.7.6
    END_VERSIONS
    """

}
