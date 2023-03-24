process FUNMAP_VIZ {
    tag  "${if (workflow.stubRun) 'funmap_viz_stub' else 'funmap_viz'}"
    label 'process_medium'
    executor 'local' //visualization will be done on local machine, connect to cytscape server running on local machine
    container 'registry.gitlab.com/bzhanglab/py4cytoscape:1.6.0'
    containerOptions '--net host' //allow the container ccess service on local host

    input:
    path funmap_el
    path ice_results
    path ice_enrich_results
    path funmap_config
    path funmap_input_data

    output:
    path 'plot_*.tsv', emit: tsv_results
    path 'plot_*.pdf', emit: pdf_results
    path 'plot_*.cys', emit: cys_results
    path 'ice_corum_overlap', emit: ice_corum_overlap
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Download input files
    wget "https://drive.google.com/uc?id=1eI4Rl167-VNv7KjsUMmjcP--lSxCoAp2&export=download" -O input.tar.gz
    tar -xzf input.tar.gz
    rm -rf input.tar.gz
    funmap_viz.py -s ${funmap_el} -i ${ice_results} -e ${ice_enrich_results} -c ${funmap_config} -d ${funmap_input_data}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       cytoscape: 3.9.1
    END_VERSIONS
    """

    stub:
    """
    wget "https://drive.google.com/uc?id=1r-FkW9crmu3_i1U7OF-MjUaDSxo9IZ5U&export=download&confirm=9iBg" -O vis_results.tgz
    tar -xzf vis_results.tgz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       cytoscape: 3.9.1
    END_VERSIONS
    """

}
