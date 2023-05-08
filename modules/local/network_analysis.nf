process NETWORK_ANALYSIS {
    tag  "${if (workflow.stubRun) 'network_analysis_stub' else 'network_analysis'}"
    label 'process_high'

    container 'registry.gitlab.com/bzhanglab/netsam:1.39.1'
    // the option '--privileged' is needed for the container to run properly
    containerOptions '--privileged'

    input:
    path funmap_el

    output:
    path 'out/*', emit: network_analysis_out
    path 'out/*.nsm', emit: netsam_nsm
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mv $funmap_el funmap.net
    export ALLOW_WGCNA_THREADS=${task.cpus}
    network_analysis.R --nthreads ${task.cpus} --input funmap.net --min-module-size ${params.min_module_size} --enrich-top-n ${params.module_enrich_top_n}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       NetSAM: 1.39.1
    END_VERSIONS
    """

    stub:
    """
    wget "https://drive.google.com/uc?id=1oFv2_SNdUoW383gamrdkdwVcDQBHQhah&export=download&confirm=9iBg" -O network_analysis_results.tgz
    tar -xzf network_analysis_results.tgz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      NetSAM: 1.39.1
    END_VERSIONS
    """
}
