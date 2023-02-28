process ICE {
    tag 'ice'
    label 'process_single'

    container 'registry.gitlab.com/bzhanglab/ice:0.1.2'

    input:
    path funmap_el

    output:
    path '*.tsv', emit: ice_results
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script: // This script is bundled with the pipeline, in mo-core/ifunmap/bin/
    """
    ice ${funmap_el} ${params.ice_min_sz} > ice_funmap_${params.ice_min_sz}.tsv
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       ice: 0.1.2
    END_VERSIONS
    """
}
