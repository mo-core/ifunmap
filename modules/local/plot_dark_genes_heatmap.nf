process PLOT_DARK_GENES_HEATMAP {
    tag  'plot_dark_genes_heatmap'
    label 'process_medium'
    container 'registry.gitlab.com/bzhanglab/python:3.8.13'
    containerOptions '--privileged'

    input:
    path funmap_el
    path dark_gene_tgi
    path gene_pubmed_count

    output:
    path '*.pdf', emit: pdf
    path '*.tsv', emit: tsv
    path 'versions.yml', emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mv ${gene_pubmed_count} gene_pubmed_count.tsv
    plot_dark_genes_heatmap.py -e ${funmap_el} -t ${dark_gene_tgi} -c gene_pubmed_count.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
       python: 3.8.13
    END_VERSIONS
    """

}
