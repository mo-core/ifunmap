#!/usr/bin/env python
import os
from os import path as osp
import argparse
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import colors
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType fonts for ease of editing in illustrator
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['savefig.dpi'] = 500

def plot_pie(df, n_dark_genes):
    ontology = df['ontology'].unique().tolist()
    ontology.sort()
    assert ontology == ['BP', 'CC', 'MF']
    ontology_full = ['Biological process', 'Molecular function',
                        'Cellular component']
    fig = plt.figure(figsize=(4, 8))
    gs0 = gridspec.GridSpec(4, 1, height_ratios=[4, 4, 4, 1], figure=fig)

    for i in range(3):
        ax = fig.add_subplot(gs0[i, 0])
        cur_ontology = ontology[i]
        cur_df = df[df['ontology'] == cur_ontology]

        # pvalue > 0.05
        # pvalue <= 0.05 & fdr > 0.05
        # fdr <= 0.05
        val = [ len(cur_df[cur_df['pvalue'] > 0.05]),
                len(cur_df[(cur_df['pvalue'] <= 0.05) &
                    (cur_df['fdr'] > 0.05)]),
            len(cur_df[cur_df['fdr'] <= 0.05])
        ]

        total = sum(val)
        # if total is less than n_dark_genes, add the difference to the first
        if total < n_dark_genes:
            val[0] += n_dark_genes - total

        def func(pct, allvalues):
            absolute = round(pct / 100.*np.sum(allvalues))
            if pct == 0:
                return ''
            return "{:.1f}%\n({:d})".format(pct, absolute)

        labels = ['pvalue > 0.05',
                'pvalue <= 0.05 and FDR > 0.05',
                'FDR <= 0.05']
        wedges, texts, autotexts = ax.pie(val, autopct =
                                        lambda pct: func(pct, val),
                                        # colors=['#fdae5b', '#8c2d04'])
                                        colors=['#fee6ce', '#fdae6b', '#8c2d04'])
        # change position of autotexts,  this is a hack
        plt.setp(autotexts, size = 16, color = 'w')
        if i == 0:
            txt = autotexts[0]
            if val[0] == 0:
                txt.set_visible(False)
            else:
                txt.set_position((1.5, 0.1))
                plt.setp(txt, size = 16, color = 'black')
            txt = autotexts[1]
            txt.set_position((0.5, 0.4))
            plt.setp(txt, size = 16, color = 'black')
            txt = autotexts[2]
            txt.set_position((-0.1, -0.4))
            plt.setp(txt, size = 16, color = 'white')
        if i == 1:
            txt = autotexts[0]
            if val[0] == 0:
                txt.set_visible(False)
            else:
                txt.set_position((1.5, 0.1))
                plt.setp(txt, size = 16, color = 'black')
            txt = autotexts[1]
            txt.set_position((0.5, 0.4))
            plt.setp(txt, size = 16, color = 'black')
            txt = autotexts[2]
            txt.set_position((-0.1, -0.4))
            plt.setp(txt, size = 16, color = 'white')
        if i == 2:
            txt = autotexts[0]
            if val[0] == 0:
                txt.set_visible(False)
            else:
                txt.set_position((1.5, 0.1))
                plt.setp(txt, size = 16, color = 'black')
            txt = autotexts[1]
            txt.set_position((0.5, 0.4))
            plt.setp(txt, size = 16, color = 'black')
            txt = autotexts[2]
            txt.set_position((-0.1, -0.4))
            plt.setp(txt, size = 16, color = 'white')


        ax.set_title(ontology_full[i], fontsize=16, fontweight='bold')
        ax.set(aspect="equal")

    # legend
    ax = fig.add_subplot(gs0[3, 0])
    ax.set_axis_off()
    ax.legend(wedges, labels, loc ="center")
    plt.savefig('dark_gene_enrich_pie.pdf', bbox_inches = 'tight', pad_inches = 0.1)


def arg_parse():
    parser = argparse.ArgumentParser(description='command line arguments.')
    parser.add_argument('-i', '--input', type=str, required=True,
                        help='path to the enrichment results file')
    parser.add_argument('-n', '--neighbor', type=str, required=True,
                        help='path to the top neighbor list')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()
    df = pd.read_csv(args.input, sep='\t')
    # read tsv without header
    top_neighbor = pd.read_csv(args.neighbor, sep='\t', header=None)
    n_dark_genes = top_neighbor.shape[1]

    all_df_g = df[['dark_gene', 'ontology', 'pValue', 'FDR']]
    all_df_g = all_df_g.groupby(['dark_gene', 'ontology']).agg(
                    pvalue=('pValue', min),
                    fdr=('FDR', min)
                )
    all_df_g = all_df_g.reset_index()
    plot_pie(all_df_g, n_dark_genes)
