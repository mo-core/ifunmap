process INPUT_FILE_CHECK {
    tag "$config_file"
    label 'process_single'

    container 'jupyter/datascience-notebook:python-3.8'

    input:
    path config_file
    path data_file

    output:
    path 'config.yml' , emit: config_file
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in mo-core/ifunmap/bin/
    """
    check_input.py \\
        -c $config_file \\
        -d $data_file

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
