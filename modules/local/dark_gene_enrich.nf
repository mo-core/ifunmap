process DARK_GENE_ENRICH {
    tag  "${if (workflow.stubRun) 'dark_gene_enrich_stub' else 'dark_gene_enrich'}"
    label 'process_high'
    container 'registry.gitlab.com/bzhanglab/webgestaltr:0.4.5'

    input:
    path funmap_el
    path gene_pubmed_count

    output:
    path 'top_neighbor.tsv', emit: top_neighbors
    path 'enrich_results.tsv', emit: enrich_results
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # first find out how many dark genes are there
    mv ${gene_pubmed_count} gene_pubmed_count.tsv
    # for now we use 4 cpus due to memory limit
    dark_gene_enrichment.R  \
        --network-el ${funmap_el} \
        --gene-pubmed-count gene_pubmed_count.tsv \
        --n-cpus 4
        # --n-cpus ${task.cpus}
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       webgestaltR: 0.4.5
    END_VERSIONS
    """

    // stub:
    // """
    // wget "https://drive.google.com/uc?id=1M92JmXsZyYqbZnxwTsn2onflbKW_msyI&export=download" -O ice_funmap_5_enrich.tsv
    // cat <<-END_VERSIONS > versions.yml
    // "${task.process}":
    //    webgestaltR: 0.4.5
    // END_VERSIONS
    // """

}
