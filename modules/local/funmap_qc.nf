process FUNMAP_QC {
    tag  "funmap_qc"
    label  "process_low"

    container 'registry.gitlab.com/bzhanglab/funmap:latest'

    input:
    path config_file
    path data_file

    output:
    path 'results/config.json', emit: funmap_config
    path 'results/figures/*', emit: funmap_figures
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    funmap qc -c ${config_file} -d ${data_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       funmap: \$(fumap --version)
    END_VERSIONS
    """

}
