process CLIQUE_ENRICH {
    //tag  "${if (workflow.stubRun) 'clique_enrich_stub' else 'clique_enrich'}"
    tag  'clique_enrich'
    label 'process_medium'

    container 'registry.gitlab.com/bzhanglab/webgestaltr:0.4.5'

    input:
    path funmap_el
    path clique_tsv

    output:
    path '*_enrich.tsv', emit: clique_enrich_results
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    wget "https://drive.google.com/uc?id=1uoGhuDlfIkloJNpeJ0sUrp2_RSBLvVx0&export=download" -O input.tar.gz
    tar -xzf input.tar.gz
    clique_complex_enrich.R  \
        --clique-file ${clique_tsv} \
        --mapping-file ./BCM_ensembl_symbol_mapping.tsv \
        --corum-gmt ./CORUM.gmt \
        --bioplex-gmt ./BioPlex3.gmt \
        --network-file ${funmap_el}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       webgestaltR: 0.4.5
    END_VERSIONS
    """
}
