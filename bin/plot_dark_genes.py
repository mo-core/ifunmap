#!/usr/bin/env python
import os
from os import path as osp
import glob
import re
import argparse
import pandas as pd
import numpy as np
import matplotlib
from matplotlib import colors
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import networkx as nx
import py4cytoscape as p4c
import seaborn as sns
from statannotations.Annotator import Annotator


matplotlib.rcParams['pdf.fonttype'] = 42  # TrueType fonts for ease of editing in illustrator
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['font.size'] = 12
matplotlib.rcParams['savefig.dpi'] = 500

def plot_heatmap(cnt_df, all_edges_dict, tgi_df, tgi_col_type):
    output_dir = '.'
    bar_color = '#1f77b4'
    n_track = len(tgi_col_type) + 1  # +1 for the edge count
    fig = plt.figure(figsize=(8, 8*5.6/(5.5+0.5*n_track)))
    gs0 = gridspec.GridSpec(3, 3, width_ratios=[1.5, 4, 0.5*n_track],
                            height_ratios=[1.5, 4, 0.1], figure=fig)

    ax00 = fig.add_subplot(gs0[0, 0])
    ax00.spines['top'].set_visible(False)
    ax00.spines['bottom'].set_visible(False)
    ax00.spines['left'].set_visible(False)
    ax00.spines['right'].set_visible(False)
    ax00.set_axis_off()
    ax00.text(0.5, 0.5, 'GeneRIF\nCount', fontsize=12, ha='center', va='center')

    ax01 = fig.add_subplot(gs0[0, 1])
    ax01.set_yscale('log')
    ax01.spines['top'].set_visible(False)
    ax01.spines['left'].set_visible(False)
    ax01.spines['right'].set_visible(False)
    # ax01.spines['bottom'].set_visible(False)
    ax01.margins(x=0)
    ax01.set_xticks([], [])
    ax01.bar(cnt_df['gene_symbol'], cnt_df['count'], color=bar_color, alpha=0.5)
    ax01.set_yticks([10, 100, 1000])
    ax01.set_yticklabels(['10', '100', '1000'], fontsize=10)

    ax10 = fig.add_subplot(gs0[1, 0])
    ax10.barh(cnt_df['gene_symbol'], np.flip(cnt_df['count'].values),
                color=bar_color, alpha=0.5)
    ax10.set_xscale('log')
    ax10.set_xticks([10, 100, 1000])
    ax10.set_xticklabels(['10', '100', '1000'], fontsize=10, rotation=90)
    ax10.invert_xaxis()
    ax10.xaxis.set_ticks_position('top')
    ax10.tick_params(axis='y', top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)
    ax10.spines['top'].set_visible(False)
    ax10.spines['left'].set_visible(False)
    # ax10.spines['right'].set_visible(False)
    ax10.spines['bottom'].set_visible(False)
    ax10.margins(y=0)
    ordered_genes = cnt_df['gene_symbol'].values
    cnt_df.to_csv(osp.join(output_dir, 'count.tsv'), index=False, sep='\t')

    # heatmap
    hist2d = True
    ax11 = fig.add_subplot(gs0[1, 1])
    n_genes = len(cnt_df['gene_symbol'])
    arr_df = pd.DataFrame(0, index=cnt_df['gene_symbol'],
                        columns=cnt_df['gene_symbol'])

    for i in range(n_genes):
        for j in range(i+1, n_genes):
            gene_p =[ordered_genes[i], ordered_genes[j]]
            gene_p.sort()
            gene_t = tuple(gene_p)
            if gene_t in all_edges_dict:
                arr_df.iloc[i, j] = 1

    arr = arr_df.values
    arr = arr + arr.T

    ax11.tick_params(axis='y', top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)
    ax11.tick_params(axis='x', top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)
    if not hist2d:
        ax11.imshow(arr, cmap='Greys',  interpolation='nearest')
    else:
        row, col = np.where(arr > 0)
        h = ax11.hist2d(row, arr.shape[0]- col, bins=500, density=True, cmap='Greens',
                        norm=colors.LogNorm())
        ax21 = fig.add_subplot(gs0[2, 1])
        ax21.set_title('density', loc='left', fontsize=8)
        cbar = plt.colorbar(h[3], cax=ax21, orientation='horizontal')
        cbar.ax.tick_params(labelsize=8)

    # barplot title section
    num_colors = 10
    color_list = list(mcolors.TABLEAU_COLORS.values())[:num_colors]
    hex_colors = [mcolors.to_hex(color) for color in color_list]

    gs02 = gridspec.GridSpecFromSubplotSpec(1, n_track, wspace = .1, hspace = .1,
                subplot_spec=gs0[0,2])
    cur_ax = fig.add_subplot(gs02[0, 0])
    cur_ax.text(0.5, 0, 'Edge count (log2)', rotation=90,
                fontsize=8, ha='center', va='bottom', color='#756bb1')
    cur_ax.axis('off')
    for i, track_name in enumerate(sorted(tgi_col_type.keys())):
        cur_ax = fig.add_subplot(gs02[0, i+1])
        cur_ax.text(0.5, 0, track_name, rotation=90,
                    fontsize=8, ha='center', va='bottom', color=hex_colors[i])
        cur_ax.axis('off')

    # now add tracks
    gs11 = gridspec.GridSpecFromSubplotSpec(1, n_track, wspace=.1, hspace=.1,
                subplot_spec=gs0[1,2])
    cur_ax = fig.add_subplot(gs11[0, 0])
    cur_ax.barh(np.flip(ordered_genes), np.log2(np.flip(arr.sum(axis=0))),
                color='#756bb1')
    for tick in cur_ax.xaxis.get_major_ticks():
        tick.label1.set_fontsize(6)
    cur_ax.margins(y=0)
    cur_ax.tick_params(axis='y', top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)

    for i, track_name in enumerate(sorted(tgi_col_type.keys())):
        cur_ax = fig.add_subplot(gs11[0, i+1])
        cur_df = pd.DataFrame({'gene_symbol': np.flip(ordered_genes)})
        cur_df = pd.merge(cur_df, tgi_df[track_name], how='left', on='gene_symbol')
        cur_df = cur_df.fillna(0)
        cur_df.loc[::-1].to_csv(osp.join(output_dir, f'{track_name}.tsv'), index=False, sep='\t')
        cur_col_type = tgi_col_type[track_name]
        if cur_col_type == 'BIN':
            cur_df[track_name] = cur_df[track_name].astype(int)
        elif cur_col_type == 'CON':
            cur_df[track_name] = cur_df[track_name].astype(float)
        cur_ax.barh(cur_df['gene_symbol'], cur_df[track_name], color=hex_colors[i])
        cur_ax.margins(x=0)
        cur_ax.margins(y=0)
        if cur_col_type == 'BIN':
            cur_ax.tick_params(axis='x', top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)
            cur_ax.tick_params(axis='y', top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)
        elif cur_col_type == 'CON':
            for tick in cur_ax.xaxis.get_major_ticks():
                tick.label1.set_fontsize(6)
                tick.label1.set_rotation(90)
            cur_ax.tick_params(axis='y', top=False, bottom=False, left=False, right=False,
                labelleft=False, labelbottom=False)

    plt.savefig(os.path.join(output_dir, 'dark_gene_heatmap.pdf'), bbox_inches = 'tight',
                pad_inches = 0.1)



def plot_dark_gene(network_el, dark_gene_tgi, gene_pubmed):
    df = pd.read_csv(network_el, sep = '\t', header=None)
    cols = df.columns
    all_edges = [tuple(sorted(x)) for x in zip(df.pop(cols[0]), df.pop(cols[1]))]
    all_edges_dict = dict.fromkeys(all_edges, 1)
    df = pd.read_csv(network_el, sep = '\t', header=None)
    all_genes_in_network = np.unique(df.values.flatten())

    g2p_df = pd.read_csv(gene_pubmed, sep = '\t')
    column_names = g2p_df.columns.tolist()
    column_names[0] = 'gene_symbol'
    column_names[1] = 'count'
    g2p_df.columns = column_names
    g2p_df = g2p_df.dropna()
    g2p_df = g2p_df.drop_duplicates(subset=['gene_symbol'], keep='first')

    cnt_df = pd.DataFrame(dict(gene_symbol=all_genes_in_network))
    cnt_df = pd.merge(cnt_df, g2p_df, how='left')
    cnt_df['count'] = cnt_df['count'].fillna(0)
    cnt_df = cnt_df.sort_values(by=['count', 'gene_symbol'], ascending=[False, True])

    tgi_df = None
    tgi_col_type = None
    # it is possible that the file is a dummy file
    if os.path.getsize(dark_gene_tgi) > 0:
        tgi_df = pd.read_csv(dark_gene_tgi, sep='\t')
        tgi_df.columns = tgi_df.columns.str.lower().str.strip()
        tgi_df.set_index('gene_symbol', inplace=True)
        tgi_col_type = tgi_df.iloc[0].to_dict()
        tgi_df = tgi_df.iloc[1:]
        if len(tgi_df.columns) > 10:
            tgi_df = tgi_df.iloc[:, :10]
            print('WARNING: more than 10 columns in dark gene TGI file. Only using the first 10 columns.')


    # cnt_df = pd.DataFrame(dict(gene_symbol=all_genes_in_network))
    # cnt_df = pd.merge(cnt_df, g2p_df, how='left')
    # cnt_df['count'] = cnt_df['count'].fillna(0)
    # cnt_df = cnt_df.sort_values(by=['count', 'gene_symbol'], ascending=[False, True])

    # cancer_driver_df = pd.read_csv(cancer_driver, sep = '\t')
    # cancer_driver_df = cancer_driver_df.drop_duplicates(subset=['gene_name'])
    # cancer_driver_df = cancer_driver_df[['gene_name']]
    # cancer_driver_df['yes'] = 1

    plot_heatmap(cnt_df, all_edges_dict, tgi_df, tgi_col_type )



def arg_parse():
    parser = argparse.ArgumentParser(description='command line arguments.')
    parser.add_argument('-e', '--edge-list-file', type=str, required=True,
                        help='path to the input network edge list file. (gene symbol)')
    parser.add_argument('-t', '--dark-gene-tgi', type=str, required=True,
                        help='path to tsi file for gene annotation.')
    parser.add_argument('-c', '--gene-pubmed-count', type=str, required=True,
                        help='path to gene pubmed count file.')
    args = parser.parse_args()
    return args

if __name__ == '__main__':
    args = arg_parse()
    plot_dark_gene(args.edge_list_file, args.dark_gene_tgi, args.gene_pubmed_count)
