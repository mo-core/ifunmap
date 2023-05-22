process PLOT_DARK_GENE_ENRICH_PIE {
    tag  'plot_dark_genes_enrich_pie'
    label 'process_medium'
    container 'registry.gitlab.com/bzhanglab/python:3.8.13'
    containerOptions '--privileged'

    input:
    path enrich_results
    path top_neighbors

    output:
    path '*.pdf', emit: pdf
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    plot_dark_gene_enrich_pie.py -i ${enrich_results} -n ${top_neighbors}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: 3.8.13
    END_VERSIONS
    """

}
