process FUNMAP {
    tag "$config_file"
    label 'process_single'

    container 'registry.gitlab.com/bzhanglab/funmap:latest'

    input:
    path config_file

    output:
    path 'config.yml' , emit: yml
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in mo-core/ifunmap/bin/
    """








    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
    END_VERSIONS
    """
}
