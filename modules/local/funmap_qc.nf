process FUNMAP_QC {
    tag  "funmap_qc"
    label  "process_low"

    container 'registry.gitlab.com/bzhanglab/funmap:latest'

    input:
    path config_file
    path data_path_file
    path data_file

    output:
    path 'results/config.json', emit: funmap_config
    path 'results/figures/*', emit: funmap_figures
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    read -r target_file < ${data_path_file}
    if [ -e "\$target_file" ]; then
        target_directory=\$(dirname "\$target_file")
        mkdir -p "\$target_directory"
        mv ${data_file} "\$target_file"

    funmap qc -c ${config_file}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       funmap: \$(fumap --version)
    END_VERSIONS
    """

   stub:
    """
    wget "https://drive.google.com/uc?id=1y6D26Hlf-Cigdz9tzrMZ3EWR55HclV9z&export=download" -O funmap_qc_results.tar.gz
    tar -xzf funmap_qc_results.tar.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       funmap: \$(funmap --version)
    END_VERSIONS
    """

}
