#!/usr/bin/env python
import argparse
import os
import pandas as pd
import csv
import warnings
from sklearn.exceptions import ConvergenceWarning
import networkx as nx
import py4cytoscape as p4c


warnings.simplefilter("ignore", category=ConvergenceWarning)

def arg_parse():
    parser = argparse.ArgumentParser(description='command line arguments.')
    parser.add_argument('-m', '--module-info', type=str, required=True,
                        help='path to the module info file')
    parser.add_argument('-l', '--module-list', type=str, required=True,
                        help='path to the module list file')
    parser.add_argument('-e', '--network-edge-list', type=str, required=True,
                        help='path to funmap edge list')
    return parser.parse_args()


def plot_modules(result_dir, network_file, module_list, module_dict):
    print(module_list)
    for mod in module_list:
        gene_list = module_dict[mod]
        funmap = nx.read_edgelist(network_file)
        subnet = funmap.subgraph(gene_list)
        edgelist_file = os.path.join(result_dir, f'{mod}_edgelist.txt')
        nodelist_file = os.path.join(result_dir, f'{mod}_genelist.txt')
        nx.write_edgelist(subnet, edgelist_file, delimiter='\t', data=False)
        subnet_nodes = pd.DataFrame.from_dict(subnet.degree)
        subnet_nodes = subnet_nodes.rename(columns={0:'node', 1:'degree'})
        subnet_nodes.to_csv(nodelist_file, sep='\t', index=False)

        edgelist_df = pd.read_csv(edgelist_file, header=None, sep='\t')
        edgelist_df = edgelist_df.rename(columns={0:'source', 1:'target'})
        node_df = pd.read_csv(nodelist_file, sep='\t')
        node_df = node_df.rename(columns={'node': 'id'})
        max_size = max(node_df['degree'])
        min_size = min(node_df['degree'])
        print(max_size, min_size)
        g = p4c.create_network_from_data_frames(node_df, edgelist_df,
                                                title=f'{mod}',
                                                collection='cytoscape_network')
        defaults = {'NODE_SHAPE': 'ellipse',
                    'EDGE_WIDTH': 1, 'EDGE_STROKE_UNSELECTED_PAINT': '#D9D9D9',
                    'NODE_FILL_COLOR': '#96CD81',
                    'NODE_BORDER_WIDTH': 0
                    }
        node_labels = p4c.map_visual_property('node label', 'id', 'p')
        node_sizes = p4c.map_visual_property('node size', 'degree', 'c',
                                            [min_size, max_size], [min_size, max_size])
        p4c.create_visual_style('my_style', defaults=defaults, mappings=[node_labels, node_sizes])
        p4c.set_visual_style('my_style')
        # input('layout manually now .... after done, press Enter to continue')
        p4c.export_image(os.path.join(result_dir, f'{mod}_cy.pdf'), type='pdf', overwrite_file=True)

    p4c.save_session(os.path.join(result_dir, f'saved_session_cy.cys'))


if __name__ == '__main__':
    args = arg_parse()
    output_dir = '.'
    module_info = args.module_info

    module_dict = {}
    with open(module_info, 'r', newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            key = row[0]
            values = row[1:]
            module_dict[key] = values

    module_list = []
    with open(args.module_list, 'r') as file:
        lines = file.readlines()
        for line in lines:
            item = line.strip()
            module_list.append(item)

    plot_modules(output_dir, args.network_edge_list, module_list, module_dict)
