process FUNMAP {
    tag  "${if (workflow.stubRun) 'funmap_stub' else 'funmap'}"
    label  "${if (workflow.stubRun) 'process_funmap_stub' else 'process_funmap'}"

    container 'registry.gitlab.com/bzhanglab/funmap:latest'

    input:
    path config_file
    path data_file

    output:
    path 'results/config.json', emit: funmap_config
    path 'results/networks/funmap.tsv', emit: funmap_el
    path 'results/networks/network_*.tsv', emit: funmap_networks
    path 'results/llr_*', emit: funmap_llr
    path 'results/figures/*', emit: funmap_figures
    path 'results/saved_*', emit: funmap_saved
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    funmap -c ${config_file} -d ${data_file} -o results

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       funmap: \$(fumap --version)
    END_VERSIONS
    """

    stub:
    """
    wget "https://drive.google.com/uc?id=1oComOLXbB_OSd4OoAI1jhoX3tqUEq-Cs&export=download" -O results.tar.gz
    tar -xzf results.tar.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       funmap: \$(funmap --version)
    END_VERSIONS
    """
}
