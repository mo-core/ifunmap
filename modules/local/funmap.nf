process FUNMAP {
    tag "funmap"
    label 'process_funmap'

    container 'registry.gitlab.com/bzhanglab/funmap:latest'

    input:
    path config_file
    path data_file

    output:
    path 'results', emit: funmap_results
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in mo-core/ifunmap/bin/
    """
    funmap -c ${config_file} -d ${data_file} -o results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       funmap: \$(fumap --version)
    END_VERSIONS
    """
}
