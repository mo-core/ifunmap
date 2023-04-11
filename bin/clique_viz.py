#!/usr/bin/env python
import os
import io
import re
import itertools
import argparse
import warnings
import random
import gc
import yaml
import tarfile
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import ConnectionPatch
from bioservices import BioMart
import igraph
import py4cytoscape as p4c
from matplotlib_venn import venn3, venn3_circles, venn3_unweighted
import networkx as nx
import seaborn as sns
import pypdf


mpl.rcParams['pdf.fonttype'] = 42  # TrueType fonts for ease of editing in illustrator
mpl.rcParams['font.family'] = 'Arial'
mpl.rcParams['font.size'] = 24
mpl.rcParams['savefig.dpi'] = 600
mpl.use('Agg')


def merge_and_delete(fig_dir, file_list, output_file):
    """
    Merge multiple PDF files into one and delete the original files.
    """
    pdf_writer = pypdf.PdfWriter()

    total_page_num = 0
    for file in file_list:
        pdf_reader = pypdf.PdfReader(os.path.join(fig_dir, file))
        cur_page_num = len(pdf_reader.pages)
        for page in range(cur_page_num):
            pdf_writer.add_page(pdf_reader.pages[page])
        pdf_writer.add_outline_item(os.path.splitext(file)[0], total_page_num)
        total_page_num = total_page_num + cur_page_num

    with open(os.path.join(fig_dir, output_file), 'wb') as fh:
        pdf_writer.write(fh)
        print('figures have been merged.')

    for filename in file_list:
        try:
            os.remove(os.path.join(fig_dir, filename))
        except:
            print(f'{filename} could not be deleted.')


def make_valid_filename(filename):
    # remove any invalid characters from the filename
    filename = re.sub(r'[^\w\s-]', '', filename)
    # replace any whitespace characters with an underscore
    filename = re.sub(r'\s+', '_', filename)
    # remove any leading or trailing hyphens
    filename = re.sub(r'^-+|-+$', '', filename)
    # return the cleaned-up filename
    return filename

def add_commas(num):
    """
    Takes a number as input and returns it as a string with commas separating the thousands.
    """
    num_str = str(num)
    if len(num_str) <= 3:
        return num_str
    else:
        return add_commas(num_str[:-3]) + ',' + num_str[-3:]


def random_colors(n=10):
    """ return a list of hex codes of lighter colors"""
    colors = []
    for i in range(n):
        rand = lambda: random.randint(127, 255)
        colors.append('#%02X%02X%02X' % (rand(), rand(), rand()))
    return colors


# https://tinyurl.com/4uv3x8x4
# there is no source code from the above link for this function
def load_corum_complexes(corum_complex_file, id_mapping=None):
    """
    Loads complexes from a CORUM complexes file and returns a pandas DataFrame.

    Parameters
    ----------
    corum_complex_file : str
        The path to the CORUM complexes file.
    id_mapping : pd.DataFrame, optional
        A pandas DataFrame containing gene symbol to Ensembl gene ID mapping.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the following columns:
            - complex_id: The CORUM complex ID.
            - complex_name: The name of the complex.
            - subunits: A set of gene symbols for the subunits of the complex. If id_mapping is not None, the gene symbols
                    are mapped to Ensembl gene IDs using the mapping.

    Examples
    --------
    >>> complexes = load_corum_complexes('corum_complexes.tsv')
    >>> print(complexes.head())

        complex_id  complex_name                       subunits
    0        10002  Cytoplasmic  {C14ORF166, CHCHD2, HAX1, HSPD1}
    1        10005  Cytoplasmic                 {ALDOA, ENO1, GAPDH}
    2        10008    Membranal      {ADAM17, CDH1, CTNNB1, GJA1}
    3        10009    Membranal     {APLP2, PSEN1, PSENEN, PSEN2}
    4        10011    Membranal       {ATP1A1, ATP1B1, ATP2B1, FXYD1}
    """
    complexes = pd.read_csv(corum_complex_file, sep='\t')
    complexes = complexes.loc[complexes['Organism'] == 'Human',
                        ['ComplexID', 'ComplexName', 'subunits(Gene name)']]
    complexes = complexes.rename(columns={
        'ComplexID': 'complex_id',
        'ComplexName': 'complex_name',
        'subunits(Gene name)': 'subunits'})
    complexes['subunits'] = complexes['subunits'].apply(lambda x: x.upper())
    if id_mapping is not None:
        complexes['subunits'] = complexes['subunits'].apply(
            lambda x: set(id_mapping.loc[id_mapping['hgnc_symbol'].isin(x.split(';')),
                                    'ensembl_gene_id'].to_list()))
    else:
        complexes['subunits'] = complexes['subunits'].apply(lambda x: set(x.split(';')))
    return complexes


def map_nw_ids(df, id_in, id_out, via=None, suffixes=('_a', '_b'),
            directed=False, space_iii=True, agg=None):
    """Tables mapping between different protein/ORF/gene identifiers.
    All combinations of pairs of IDs are returned. Mappings are symmetric.
    Supported IDs are: orf_id, ensembl_gene_id, uniprot_ac
    Where orf_id refers to our internal ORF IDs.
    The ORF ID mappings come from the ORFeome annotation project where ORFs
    were sequence aligned to ensembl transcripts and proteins. The uniprot
    to ensembl mappings are from a file provided by UniProt.
    Note:
        The ensembl gene IDs and UniProt ACs are an unfiltered redundant set.
        You will probably want to filter the set after mapping.
    Args:
        df (DataFrame): table of protein-protein interactions,
                        one row for each pair
        id_in/id_out (str): gene/protein identifiers to map between
        via (str): optional identifier to use as an intermediate in the
                mapping. Useful for example when mapping uniprot IDs to
                our ORF IDs from experiments that can ony determine the
                gene-level then can map via ensembl_gene_id.
        suffixes (tuple(str, str)): of the two protein identifiers of each pair
        directed (bool): if True, retain the mapping of the two proteins,
                        if False, sort the IDs and drop duplicates
        space_iii (bool): restrict ORF IDs to Space III, where there is one
                            ORF per gene. Only has an effect if one of id_in,
                            id_out or via_ID is 'orf_id'.
        agg (function): Optional. Custom aggregation function to choose a
                        single pair or combine multiple pairs when input pairs
                        are mapped to the same pair in the new ID. Function
                        must take a DataFrame and return a DataFrame with no
                        duplicate pairs.
    Returns:
        DataFrame: PPI dataset mapped to new protein/gene ID. All unique
                combinations of the new ID are mapped to, so one pair in
                the input dataset can map to multiple pairs in the output
                and vice-versa.
    """
    id_in_a = id_in + suffixes[0]
    id_in_b = id_in + suffixes[1]
    id_out_a = id_out + suffixes[0]
    id_out_b = id_out + suffixes[1]
    id_map = load_id_map(id_in, id_out)
    out = df.copy()
    out = (pd.merge(out, id_map,
                    how='inner',
                    left_on=id_in_a,
                    right_on=id_in)
                .drop(columns=id_in)
                .rename(columns={id_out: id_out_a}))
    out = (pd.merge(out, id_map,
                    how='inner',
                    left_on=id_in_b,
                    right_on=id_in)
                .drop(columns=id_in)
                .rename(columns={id_out: id_out_b}))
    if out.loc[:, [id_out_a, id_out_b]].isnull().any().any():
        raise UserWarning('Unexpected missing values')
    if not directed:
        a = out[[id_out_a, id_out_b]].min(axis=1)
        b = out[[id_out_a, id_out_b]].max(axis=1)
        out[id_out_a] = a
        out[id_out_b] = b
    out = out.drop(columns=[id_in_a, id_in_b])
    pair_duplicates = out.duplicated(subset=[id_out_a, id_out_b], keep=False)
    if (out.duplicated(keep=False) != pair_duplicates).any() and agg is None:
        warnings.warn('Warning: mapping between gene/protein identifiers has '
                        'resulted in different pairs in the input ID being mapped to '
                        'the same pair in the output ID.\n'
                        'You may wish to use the `agg` argument to customize '
                        'the choice of which of the pair\'s infomation to keep or how '
                        'to combine the information from multiple pairs.')
    if agg is None:
        out = out.drop_duplicates(subset=[id_out_a, id_out_b])
    else:
        out = agg(out)
        if out.duplicated().any():
            raise ValueError('Problem with your agg function')
    cols = out.columns.tolist()
    cols = cols[-2:] + cols[:-2]
    out = out.loc[:, cols]
    out = out.set_index(out[id_out_a].astype(str) +
                        '_' +
                        out[id_out_b].astype(str))
    return out


def load_protein_coding_genome():
    return set(pd.read_csv('HuRI_nature_paper_2020_Supplementary_Table_1.txt', sep='\t')
                ['Ensembl_gene_id'].values)


def _guess_id_type_and_suffixes(columns):
    id_a, id_b = columns[:2]
    if id_a.split('_')[:-1] == id_b.split('_')[:-1]:
        return ('_'.join(id_a.split('_')[:-1]),
                ('_' + id_a.split('_')[-1],
                '_' + id_b.split('_')[-1]))
    else:
        raise UserWarning('Failed to guess id_type and suffixes. Please specify.')


def _guess_node_id(columns, suffixes=('_a', '_b')):
    """
    Given a table of edges, try and guess the node IDs.
    Should I give a warning? i.e. guessing gene id type is x, to silence pass id_type=x?
    Maybe only give warning if it's ambiguous? i.e. there are not multiple?
    """
    cols_a = [c[:-len(suffixes[0])] for c in columns if c.endswith(suffixes[0])]
    cols_b = [c[:-len(suffixes[1])] for c in columns if c.endswith(suffixes[1])]
    for col_a in cols_a:
        for col_b in cols_b:
            if col_a == col_b:
                return col_a
    raise UserWarning('Could not guess node IDs from: ' + ' | '.join(columns))


def split_node_and_edge_tables(edges, id_type=None, suffixes=None):
    es = edges.copy()
    if id_type is None:
        if suffixes is None:
            id_type, suffixes = _guess_id_type_and_suffixes(es.columns)
        else:
            id_type = _guess_node_id(es.columns, suffixes)
    elif suffixes is None:
        suffixes = es.columns[0].replace(id_type, ''), es.columns[1].replace(id_type, '')
    columns_a = [c for c in es.columns if c.endswith(suffixes[0])]
    columns_b = [c for c in es.columns if c.endswith(suffixes[1])]
    ns_a = (es.loc[:, columns_a]
            .rename(columns={c: c[:-len(suffixes[0])] for c in columns_a})
            .drop_duplicates())
    ns_b = (es.loc[:, columns_b]
            .rename(columns={c: c[:-len(suffixes[1])] for c in columns_b})
            .drop_duplicates())
    ns = (pd.concat([ns_a, ns_b])
            .drop_duplicates()
            .set_index(id_type, verify_integrity=True))
    gene_id_columns = [id_type + s for s in suffixes]
    es = es.loc[:, gene_id_columns + [c for c in es.columns if c not in (columns_a + columns_b)]]
    return ns, es


def merge_node_and_edge_tables(nodes,
                                edges,
                                id_type=None,
                                suffixes=('_a', '_b'),
                                node_id_column=None):
    """Combine data on nodes and edges into a table of edges.
    Args:
        nodes (pandas.DataFrame): table of nodes
        edges (pandas.DataFrame): table of edges
    """
    if id_type is None:
        _guess_node_id(edges.columns, suffixes)
    if node_id_column is None:
        df = pd.merge(edges,
                        nodes,
                        right_index=True,
                        left_on=id_type + suffixes[0],
                        how='left')
        df = pd.merge(df,
                        nodes,
                        right_index=True,
                        left_on=id_type + suffixes[1],
                        how='left',
                        suffixes=suffixes)
    else:
        df = pd.merge(edges,
                        nodes,
                        right_on=node_id_column,
                        left_on=id_type + suffixes[0],
                        how='left')
        df = pd.merge(df,
                        nodes,
                        right_on=node_id_column,
                        left_on=id_type + suffixes[1],
                        how='left',
                        suffixes=suffixes)
    if df.columns.duplicated().any():
        df = df.loc[:, ~df.columns.duplicated()]
    return df


def format_network(data, fmt, id_type=None, suffixes=('_a', '_b'), directed=False):
    """Convert data format of network between pandas/networkx/igraph etc.
    fmt options are:
        - `pandas`: DataFrame
        - `igraph`: igraph Graph
        - `list`: list of a dict for each PPI
    Args:
        data (pandas.DataFrame / igraph.Graph / list): input network
        fmt (str): format to convert to; pandas/nx/igraph/list
        id_type (str): gene/protein identifier type used
        suffixes (tuple(str, str)): at the end of id_type to distinguish the two nodes
        directed (bool): return directed for igraph/networkx
    Returns:
        (pandas.DataFrame / igraph.Graph / list): network in specified format
    """
    valid_fmt = ['pandas', 'igraph', 'list']
    if fmt not in valid_fmt:
        raise UserWarning('Unsupported fmt: ' + fmt +
                        '\nValid options are: ' + '/'.join(valid_fmt))
    fmt_out = fmt
    fmt_in = 'unknown'
    for fmt_desc, data_type in [('pandas', pd.DataFrame),
                                ('igraph', igraph.Graph),
                                ('list', list)]:
        if isinstance(data, data_type):
            fmt_in = fmt_desc
            break
    if fmt_in == 'unknown':
        raise ValueError('Unsupported input type: ' + type(data))
    if id_type is None:
        if fmt_in == 'pandas':
            id_type, suffixes = _guess_id_type_and_suffixes(data.columns)
        else:
            raise ValueError('Require value for id_type argument')
    id_a = id_type + suffixes[0]
    id_b = id_type + suffixes[1]

    if fmt_in != 'pandas' and fmt_out != 'pandas':
        # via pandas.DataFrame, so don't have to code every possible conversion
        tbl = format_network(data,
                            'pandas',
                            id_type=id_type,
                            suffixes=suffixes,
                            directed=directed)
        return format_network(tbl,
                            fmt_out,
                            id_type=id_type,
                            suffixes=suffixes,
                            directed=directed)
    elif fmt_in == 'pandas' and fmt_out == 'pandas':
        return data
    elif fmt_in == 'pandas' and fmt_out == 'igraph':
        node_df, edge_df = split_node_and_edge_tables(data, id_type=id_type, suffixes=suffixes)
        g = igraph.Graph()
        g = g.TupleList([(a, b) for a, b in edge_df[[id_a, id_b]].values],
                        directed=directed)
        for column in edge_df.columns:
            if column not in [id_a, id_b]:
                g.es[column] = edge_df[column].values
        for column in node_df.columns:
            g.vs[column] = node_df.loc[g.vs['name'], column].values
        return g
    elif fmt_in == 'pandas' and fmt_out == 'list':
        d = data.to_dict()
        mylist = [{k: d[k][idx] for k in d.keys()} for idx in d[id_a].keys()]
        return mylist
    elif fmt_in == 'igraph' and fmt_out == 'pandas':
        d = {id_a: [data.vs[e.source]['name'] for e in data.es],
            id_b: [data.vs[e.target]['name'] for e in data.es]}
        d.update({k: data.es[k] for k in data.es.attribute_names()})
        edge_df = pd.DataFrame(d)
        edge_df = edge_df.set_index(edge_df[id_a] + '_' + edge_df[id_b])
        node_df = pd.DataFrame({k: data.vs[k] for k in data.vs.attribute_names()}).set_index('name')
        out = merge_node_and_edge_tables(node_df, edge_df, id_type=id_type, suffixes=suffixes)
        return out
    elif fmt_in == 'list' and fmt_out == 'pandas':
        out = pd.DataFrame(data)
        out = out.set_index(out[id_a] + '_' + out[id_b])
        return out
    else:
        raise UserWarning('Something went wrong...')


def load_nw_bioplex(id_type='ensembl_gene_id',
                    out_file='BioPlex_3_edge_list_gene_symbol.tsv',
                    fmt='pandas'):
    """
    Load BioPlex interaction data.

    Parameters
    ----------
    id_type : str, optional
        Type of identifier to use for genes. Possible values are 'ensembl_gene_id', 'orf_id', or 'gene_symbol'.
        The default value is 'ensembl_gene_id'.
    out_file : str, optional
        Path to the output file. The default value is 'BioPlex_3_edge_list_gene_symbol.tsv'.
    fmt : str, optional
        Format to use for the output. Possible values are 'pandas', 'nx_graph', or 'edge_list'.
        The default value is 'pandas'.

    Returns
    -------
    pandas.DataFrame or tuple
        BioPlex interaction data in the specified format.

    Raises
    ------
    ValueError
        If the specified identifier type is not supported.

    Notes
    -----
    BioPlex is a high-throughput protein-protein interaction dataset. This function loads the data from the file
    'BioPlex_293T_HCT116_combined_Dec_2019.tsv' and preprocesses it to remove invalid interactions and convert
    UniProt accession numbers to gene identifiers. The resulting data is returned in the specified format.

    """

    fpath = 'BioPlex_293T_HCT116_combined_Dec_2019.tsv'
    suffixes = ('_a', '_b')

    # Learnt that A is bait and B is prey in this file from email
    # communication between Katja and Ed Huttlin. Caveat is that
    # the pairs are unique in that file so in cases where both ways
    # are found (bait-X/prey-Y & bait-Y/prey-X), then only one is chosen.
    col_bait = 'UniprotA'
    col_prey = 'UniprotB'
    col_symbol_a = 'SymbolA'
    col_symbol_b = 'SymbolB'
    bp = (pd.read_csv(fpath,
                    sep='\t',
                    header=0,
                    usecols=[col_bait, col_prey, col_symbol_a, col_symbol_b])
        .rename(columns={col_bait: 'uniprot_ac_bait',
                            col_prey: 'uniprot_ac_prey'}))
    bp = bp.loc[(bp['uniprot_ac_bait'] != 'UNKNOWN') &
                (bp['uniprot_ac_prey'] != 'UNKNOWN'), :]
    # converting from isoform level UniProt AC to protein/gene level
    bp['uniprot_ac_bait'] = bp['uniprot_ac_bait'].str.replace('-.*', '', regex=True)
    bp['uniprot_ac_prey'] = bp['uniprot_ac_prey'].str.replace('-.*', '', regex=True)
    bp['uniprot_ac_a'] = (bp.loc[:, ['uniprot_ac_bait',
                                        'uniprot_ac_prey']]
                            .min(axis=1))
    bp['uniprot_ac_b'] = (bp.loc[:, ['uniprot_ac_bait',
                                        'uniprot_ac_prey']]
                            .max(axis=1))
    bp = bp.drop(['uniprot_ac_bait', 'uniprot_ac_prey'], axis=1)
    bp = bp.drop_duplicates()
    bp = bp.set_index(bp['uniprot_ac_a'] + '_' + bp['uniprot_ac_b'])
    bp = bp.drop_duplicates()
    if id_type == 'orf_id':
        bp = bp.drop(['SymbolA', 'SymbolB'], axis=1)
        bp = map_nw_ids(bp, 'uniprot_ac', 'orf_id', via='ensembl_gene_id',
                        suffixes=suffixes, directed=False)
    elif id_type == 'ensembl_gene_id':
        bp = bp.drop(['SymbolA', 'SymbolB'], axis=1)
        bp = map_nw_ids(bp, 'uniprot_ac', 'ensembl_gene_id',
                        suffixes=suffixes, directed=False)
        pcg = load_protein_coding_genome()
        bp = bp.loc[bp['ensembl_gene_id' + suffixes[0]].isin(pcg) &
                    bp['ensembl_gene_id' + suffixes[1]].isin(pcg), :]
    elif id_type == 'gene_symbol':
        bp = bp.drop(['uniprot_ac_a', 'uniprot_ac_b'], axis=1)
        bp = bp.rename(columns={'SymbolA': 'gene_symbol_a',
                                'SymbolB': 'gene_symbol_b'})
    else:
        raise ValueError('Unsupported ID type: ' + id_type)

    bp = format_network(bp,
                        fmt,
                        id_type=id_type,
                        suffixes=suffixes,
                        directed=False)

    if not os.path.exists(out_file):
        print('Writing BioPlex edge list to: ' + out_file)
        bp.to_csv(out_file, sep='\t', index=False, header=False)

    return bp


def ensembl_human_id_mappings(release_date='may2021'):
    """
    Retrieves Ensembl to external ID mappings for human genes.

    Parameters
    ----------
    release_date : str, optional
        The release date of the Ensembl database to use, in the format 'monthyear',
        e.g. 'may2021'. Default is 'may2021'.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame containing the Ensembl gene IDs, Ensembl transcript IDs,
        HGNC IDs, HGNC symbols, and Entrez gene IDs for human genes. Rows containing
        missing data are removed, and the `entrez_gene_id` column is cast to an
        `Int64` type.

    Raises
    ------
    UserWarning
        If the name of the Entrez gene ID field cannot be found in the BioMart attributes.

    """
    bm = BioMart(host=release_date + '.archive.ensembl.org')
    bm.datasets('ENSEMBL_MART_ENSEMBL')
    bm.add_dataset_to_xml('hsapiens_gene_ensembl')
    bm.add_attribute_to_xml('ensembl_gene_id')
    bm.add_attribute_to_xml('ensembl_transcript_id')
    # can have maximum 3 external identifiers
    bm.add_attribute_to_xml('hgnc_id')
    bm.add_attribute_to_xml('hgnc_symbol')
    if 'entrezgene' in bm.attributes(dataset='hsapiens_gene_ensembl').keys():
        bm.add_attribute_to_xml('entrezgene')
    elif 'entrezgene_id' in bm.attributes(dataset='hsapiens_gene_ensembl').keys():
        bm.add_attribute_to_xml('entrezgene_id')
    else:
        raise UserWarning('couldnt find correct name for entrez gene id field')
    res = bm.query(bm.get_xml())
    df = pd.read_csv(io.StringIO(res),
                        names=['ensembl_gene_id',
                            'ensembl_transcript_id',
                            'hgnc_id',
                            'hgnc_symbol',
                            'entrez_gene_id'],
                        sep='\t')
    df['entrez_gene_id'] = df['entrez_gene_id'].astype('Int64')
    return df


def load_id_map(id_in, id_out, ensembl_version='v104', release_date='may2021'):
    """
    Loads a mapping of gene IDs from one type to another.

    Parameters
    ----------
    id_in : str
        The type of gene ID to map from. Must be one of the following:
        'ensembl_gene_id', 'uniprot_ac', 'hgnc_id', 'hgnc_symbol'.
    id_out : str
        The type of gene ID to map to. Must be one of the following:
        'ensembl_gene_id', 'uniprot_ac', 'hgnc_id', 'hgnc_symbol'.
    ensembl_version : str, optional
        The version of Ensembl gene ID mappings to use. Default is 'v104'.
    release_date : str, optional
        The release date of Ensembl gene ID mappings to use. Default is 'may2021'.

    Returns
    -------
    pd.DataFrame
        A DataFrame with two columns, the input gene ID and the output gene ID.

    Raises
    ------
    ValueError
        If either id_in or id_out is not a valid gene ID type, or if id_in == id_out.

    Notes
    -----
    The supported gene ID types are 'ensembl_gene_id', 'uniprot_ac', 'hgnc_id', and 'hgnc_symbol'.
    If mapping from uniprot_ac to ensembl_gene_id, the function reads a pre-existing file.
    Otherwise, it reads Ensembl's gene ID mappings from a file or downloads them if the file does not exist.

    """
    valid_ids = ['ensembl_gene_id',
                    'uniprot_ac',
                    'hgnc_id',
                    'hgnc_symbol']
    for id_type in [id_in, id_out]:
        if id_type not in valid_ids:
            raise ValueError('Unsupported ID: ' + id_type +
                                '\nChoices are: ' + ', '.join(valid_ids))
    if id_in == id_out:
        raise ValueError('Invalid arguments: id_in == id_out')
    ids = set([id_in, id_out])
    if ids == {'uniprot_ac', 'ensembl_gene_id'}:
        return pd.read_csv('uniprot_to_ensembl.tsv', sep='\t')
    id_map_path = f'ensembl_{ensembl_version}_gene_id_mapping.tsv'
    if os.path.exists(id_map_path):
        df = pd.read_csv(id_map_path, sep='\t')
    else:
        df = ensembl_human_id_mappings(release_date=release_date)
        df.to_csv(id_map_path, sep='\t', index=False)
    df = df.loc[:, [id_in, id_out]].dropna().drop_duplicates()
    if id_in == 'hgnc_id' or id_out == 'hgnc_id':
        df['hgnc_id'] = df['hgnc_id'].str.replace('HGNC:', '').astype(int)
    if id_in == 'entrez_gene_id' or id_out == 'entrez_gene_id':
        df['entrez_gene_id'] = df['entrez_gene_id'].astype('Int64')
    return df


def arg_parse():
    parser = argparse.ArgumentParser(description='command line arguments.')
    parser.add_argument('-s', '--input-network-symbol', type=str, required=True,
                        help='path to the input network edge list file.(gene symbol)')
    parser.add_argument('-i', '--ice-module-file', type=str, required=True,
                        help='path to the output of ICE module file.')
    parser.add_argument('-e', '--enrich-results-file', type=str, required=True,
                        help='path to the enrichment results file.')
    parser.add_argument('-c', '--funmap-config-file', type=str, required=True,
                        help='path to the config file for funmap.')
    parser.add_argument('-d', '--data-file', type=str, required=True,
                        help='path to the data file for funmap.')
    args = parser.parse_args()
    return args


def plot_venn(set_data, labels, filename):
    # data is a tuple of (Abc, aBc, ABc, abC, AbC, aBC, ABC)
    # modify labels to add total sum of each set
    labels = list(labels)
    print(set_data)
    set_a_sum = set_data[0] + set_data[2] + set_data[4] + set_data[6]
    set_b_sum = set_data[1] + set_data[2] + set_data[5] + set_data[6]
    set_c_sum = set_data[3] + set_data[4] + set_data[5] + set_data[6]
    sums = [set_a_sum, set_b_sum, set_c_sum]
    print(sums)
    new_labels = []
    for (i, l) in enumerate(labels):
        new_labels.append(f'{l} ({add_commas(sums[i])})')

    plt.figure(figsize=(15,15))
    v = venn3(subsets=set_data, set_labels=new_labels)
    c = venn3_circles(subsets=set_data, linestyle='solid', linewidth=1, color='black')
    ratio_1 = (set_data[2] + set_data[6]) / set_a_sum * 100
    # round to 2 decimal places
    ratio_1 = round(ratio_1, 2)
    text_1 = f'({set_data[2]}+{set_data[6]})/{set_a_sum}={ratio_1}%'
    plt.annotate(text_1, xy=v.get_label_by_id('110').get_position() - np.array([0, 0.05]), xytext=(-250,-500),
                ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=-0.5',color='gray'))
    ratio_2 = (set_data[5] + set_data[6]) / set_c_sum * 100
    # round to 2 decimal places
    ratio_2 = round(ratio_2, 2)
    text_2 = f'({set_data[5]}+{set_data[6]})/{set_c_sum}={ratio_2}%'
    plt.annotate(text_2, xy=v.get_label_by_id('011').get_position() - np.array([-0.05, -0.02]), xytext=(200,-500),
                ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='gray', alpha=0.1),
                arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.5',color='gray'))
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def plot_ice_modules(network_edge_list, ice_module_file, base_url, style, output_file):
    p4c.networks.delete_all_networks(base_url=base_url)
    g = igraph.Graph.Read_Ncol(network_edge_list, directed=False)
    print(g.summary())
    clusters = g.connected_components(mode='weak')
    nclusters = len(clusters)
    print(f'number of conneted components: {nclusters}')
    edge_tuples = []
    clique_f = open(ice_module_file)
    gene_clique_mapping = {}
    clique_id = 0
    clique_id_for_overlap_node = 9999999
    for line in clique_f:
        genes = line.strip().split('\t')
        for gene in genes:
            if gene in gene_clique_mapping:
                gene_clique_mapping[gene] = clique_id_for_overlap_node
            else:
                gene_clique_mapping[gene] = clique_id
        for (idx, gene) in enumerate(genes):
            for j in range(idx+1, len(genes)):
                edge_tuples.append((gene, genes[j]))
        clique_id = clique_id + 1
    clique_f.close()
    # random colors
    colors = random_colors(clique_id)
    node_color_keys = []
    node_color_vals = []
    for i in range(clique_id):
        node_color_keys.append(str(i))
        node_color_vals.append(colors[i%len(colors)])
    node_color_keys.append(str(clique_id_for_overlap_node))
    node_color_vals.append('#FDF8B9')

    print('creating network')
    g = igraph.Graph.TupleList(edge_tuples, directed=False)
    for v in g.vs:
        g.vs.select(name=v['name'])['clique_id'] = gene_clique_mapping[v['name']]
    print(g.summary())
    clusters = g.connected_components(mode='weak')
    nclusters = len(clusters)
    print(f'number of conneted components: {nclusters}')
    network = p4c.networks.create_network_from_igraph(g, title='plot_funmap_ice_module', base_url=base_url)
    node_color_map = p4c.style_mappings.map_visual_property('node fill color',
                'clique_id', 'd', node_color_keys, node_color_vals, base_url=base_url)
    node_label = p4c.style_mappings.map_visual_property('node label', 'id', 'p', base_url=base_url)
    style_name = 'my_style'
    p4c.styles.create_visual_style(style_name,
                defaults=style,
                mappings=[node_color_map, node_label], base_url=base_url)
    print('laying out network')
    p4c.layouts.layout_network(layout_name='force-directed', network=network, base_url=base_url)
    print('setting visual style')
    p4c.styles.set_visual_style(style_name, network=network, base_url=base_url)
    print('exporting image ...')
    p4c.network_views.fit_content(network=network, base_url=base_url)
    p4c.network_views.export_image(filename=output_file, type='PDF', overwrite_file=True, base_url=base_url)
    print('saving session ...')
    # saving seems to be working only with absolute path
    # replace the ouput_file with .cys extension
    cys_out = output_file.replace('.pdf', '.cys')
    p4c.session.save_session(filename=os.path.abspath(cys_out), base_url=base_url)


def plot_heatmap(clique_id, funmap_nodes, max_sample, cfg_file, data_dir, out_dir):
    with open(cfg_file, 'r') as fh:
        cfg_dict = yaml.load(fh, Loader=yaml.FullLoader)

    data_files = cfg_dict['data_files']
    # only plot the one with highest pair wise correlation
    max_avg_corr = -1

    # save the corr matrix
    corr_res = {}
    for item in data_files:
        cur_name = item['name']
        cur_path = item['path']
        cur_data = pd.read_csv(os.path.join(data_dir, cur_path), sep='\t')
        selected = cur_data[cur_data[cur_data.columns[0]].isin(funmap_nodes)].copy()
        # continue if not all nodes are found
        if selected.shape[0] != len(funmap_nodes):
            continue
        print(selected.shape)
        selected.sort_values(by=selected.columns[0], inplace=True)
        selected.set_index(selected.columns[0], inplace=True)
        selected_t = selected.T
        cc = selected_t.corr().values
        n = cc.shape[0]
        cc_avg = (cc.sum() - n) / (n * (n - 1))
        corr_res[cur_name] = cc_avg
        if cc_avg > max_avg_corr:
            selected_data_file = cur_name
            max_avg_corr = cc_avg
        print(f'{cur_name}: {cc_avg}')
    gc.collect()
    corr_df = pd.DataFrame.from_dict(corr_res, orient='index', columns=['avg_corr'])
    corr_df.to_csv(os.path.join(out_dir, f'clique_{clique_id}_avg_corr.tsv'), sep='\t')
    print(f'clique {clique_id}: max avg corr: {selected_data_file}')

    plot_file_names = []
    for item in data_files:
        cur_name = item['name']
        cur_path = item['path']
        with open(os.path.join(data_dir, cur_path), 'r') as f:
            cur_data = pd.read_csv(f, sep='\t')
        selected = cur_data[cur_data[cur_data.columns[0]].isin(funmap_nodes)].copy()
        print(selected.shape)
        if selected.shape[0] != len(funmap_nodes):
            del cur_data, selected
            gc.collect()
            continue
        selected.sort_values(by=selected.columns[0], inplace=True)
        selected.set_index(selected.columns[0], inplace=True)
        # continue if not all nodes are found
        if selected.shape[1] > max_sample:
            del cur_data, selected
            gc.collect()
            continue
        selected_t = selected.T
        selected_z = ((selected_t - selected_t.mean())/selected_t.std()).T
        # sort column by mean
        selected_z_sorted = selected_z.reindex(selected_z.mean().sort_values().index, axis=1)
        fig, ax = plt.subplots(figsize=(15, 6))
        heatmap = sns.heatmap(selected_z_sorted, linewidth=0.5, cmap='RdBu_r',
                    ax=ax, vmin=-3, vmax=3)
        cbar = ax.collections[0].colorbar
        cbar.set_ticks([-3, -2, -1, 0, 1, 2, 3])
        ax.set_title(cur_name)
        ax.set_ylabel('')
        ax.tick_params(left=False, bottom=False)
        ax.xaxis.set_visible(False)
        fig.tight_layout()

        file_name = make_valid_filename(cur_name)
        cur_file_name = f'clique_{clique_id}_heatmap_{file_name}.pdf'
        fig.savefig(os.path.join('.', os.path.join(out_dir, cur_file_name)))
        plot_file_names.append(cur_file_name)
        fig.clf()
        heatmap.figure.clf()
        plt.clf()
        plt.close('all')
        del cur_data, selected, heatmap, fig, ax, cbar, selected_t, selected_z, selected_z_sorted
        gc.collect()

    return plot_file_names


def plot_ice_corum_overlap(enrich_results_file, default_style, base_url, output_dir,
                        config_file, data_file,
                        heatmap_max_size=10, heatmap_max_sample=1000):
    data_dir = './funmap_data'
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)
    with tarfile.open(data_file, "r:gz") as tar:
        for member in tar.getmembers():
            member.name = os.path.basename(member.name)
            if member.isfile():
                tar.extract(member, path=data_dir)

    p4c.networks.delete_all_networks(base_url=base_url)
    p4c.styles.delete_all_visual_styles(base_url=base_url)
    # read the enrich results
    enrich_results = pd.read_csv(enrich_results_file, sep='\t', na_filter=False)
    # iterate over the rows
    for _, row in enrich_results.iterrows():
        clique_id = row['clique_id']
        print(f'plotting clique {clique_id} ...')
        # get the clique_list column
        funmap_nodes = row['clique_list']
        # split into a set of nodes, separated by ;
        funmap_nodes = set(funmap_nodes.split(';'))
        # get the corum_geneSet column
        corum_genes = row['corum_geneSet']
        print(corum_genes)
        if corum_genes == 'NA':
            corum_nodes = set()
        else:
            # get corum_complex_list column and split into a set of nodes, separated by ;
            corum_nodes = row['corum_complex_list']
            corum_nodes = set(corum_nodes.split(';'))

        all_nodes = list(funmap_nodes | corum_nodes)
        funmap_edges = set(tuple(sorted(list(t))) for t in itertools.combinations(funmap_nodes,2))
        if len(corum_nodes) > 0:
            corum_edges = set(tuple(sorted(list(t))) for t in itertools.combinations(corum_nodes,2))
        else:
            corum_edges = set()
        all_edges = list(funmap_edges | corum_edges)
        all_nodes_type = []
        for i in all_nodes:
            if i in funmap_nodes and i in corum_nodes:
                node_type = '3'
            elif i in funmap_nodes:
                node_type = '1'
            elif i in corum_nodes:
                node_type = '2'
            all_nodes_type.append(node_type)
        node_is_in_funmap = []
        for i in all_nodes:
            if i in funmap_nodes:
                node_in_funmap = '1'
            else:
                node_in_funmap = '0'
            node_is_in_funmap.append(node_in_funmap)
        all_edges_type = []
        for i in all_edges:
            if i in funmap_edges and i in corum_edges:
                edge_type = '3'
            elif i in funmap_edges:
                edge_type = '1'
            elif i in corum_edges:
                edge_type = '2'
            all_edges_type.append(edge_type)
        edge_is_in_funmap = []
        for i in all_edges:
            if i in funmap_edges:
                edge_in_funmap = '1'
            else:
                edge_in_funmap = '0'
            edge_is_in_funmap.append(edge_in_funmap)

        g = igraph.Graph()
        g.add_vertices(all_nodes)
        g.add_edges(all_edges)
        g.vs['node_type'] = all_nodes_type
        g.vs['node_in_funmap'] = node_is_in_funmap
        g.es['edge_type'] = all_edges_type
        g.es['edge_in_funmap'] = edge_is_in_funmap
        g.vs['symbols'] = all_nodes

        network = p4c.networks.create_network_from_igraph(g, title=f'clique_{clique_id}', base_url=base_url)
        node_label = p4c.style_mappings.map_visual_property('node label', 'symbols', 'p', base_url=base_url)
        node_type = p4c.style_mappings.map_visual_property('node fill color',
                    'node_type', 'd', ['1', '2', '3'],
                    ['#e78ac3', '#00aeef', '#a29def'], base_url=base_url)
        edge_color = p4c.style_mappings.map_visual_property('edge stroke unselected paint',
                    'edge_type', 'd', ['1', '2', '3'],
                    ['#e78ac3', '#00aeef', '#a29def'], base_url=base_url)
        edge_line_style = p4c.style_mappings.map_visual_property('edge line type',
                    'edge_in_funmap', 'd', ['0', '1'],
                    ['LONG_DASH', 'SOLID'], base_url=base_url)
        style_name = 'my_style'
        p4c.styles.create_visual_style(style_name,
                    defaults=default_style,
                    mappings=[node_label, node_type, edge_color, edge_line_style], base_url=base_url)
        p4c.py4cytoscape_tuning.CATCHUP_NETWORK_TIMEOUT_SECS = 600
        p4c.layouts.layout_network(layout_name='force-directed', network=network, base_url=base_url)
        p4c.styles.set_visual_style(style_name, network=network, base_url=base_url)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        pdf_name = f'clique_{clique_id}.pdf'
        p4c.network_views.fit_content(network=network, base_url=base_url)
        p4c.network_views.export_image(filename=os.path.join(output_dir, pdf_name), type='PDF', overwrite_file=True, base_url=base_url)

        if len(funmap_nodes) <= heatmap_max_size:
            heatmap_files = plot_heatmap(clique_id, funmap_nodes, heatmap_max_sample, config_file, data_dir, output_dir)
        else:
            heatmap_files = []

        all_pdfs = [pdf_name] + heatmap_files
        merge_and_delete(output_dir, all_pdfs, f'clique_{clique_id}_all.pdf')
        os.rename(os.path.join(output_dir, f'clique_{clique_id}_all.pdf'), os.path.join(output_dir, pdf_name))

    output_cys = os.path.join(output_dir, f'ice_corum_overlap.cys')
    p4c.session.save_session(filename=os.path.abspath(output_cys), base_url=base_url)


def plot_pie_chart(results_file, output_file):
    # read tsv file
    results = pd.read_csv(results_file, sep='\t')
    # get the count of the column 'clique_class'
    counts = results['clique_class'].value_counts()
    counts = pd.DataFrame(counts)
    counts = counts.sort_index()
    # two number for the first pie chart
    pie1_1 = counts.iloc[0].values[0]
    # sum the counts of class 2, 3, 4 (index)
    pie1_2 = counts.iloc[1:].sum().values[0]

    pie2_1 = counts.iloc[1].values[0]
    pie2_2 = counts.iloc[2:].sum().values[0]

    pie3_1 = counts.iloc[2].values[0]
    pie3_2 = counts.iloc[3].values[0]

    print(pie1_1, pie1_2, pie2_1, pie2_2, pie3_1, pie3_2)

    # make figure and assign axis objects
    fig = plt.figure(figsize=(15, 4))
    ax1 = fig.add_subplot(131)
    ax2 = fig.add_subplot(132)
    ax3 = fig.add_subplot(133)
    fig.subplots_adjust(wspace=0)

    # large pie chart parameters
    ratio_11 = pie1_1 / (pie1_1 + pie1_2)
    ratio_12 = pie1_2 / (pie1_1 + pie1_2)
    ratios = [ratio_12, ratio_11]
    labels = ['', '']
    explode = [0.1, 0]
    # rotate so that first wedge is split by the x-axis
    angle = -180 * ratios[0]

    def autopct_format(values):
        def my_format(pct):
            total = sum(values)
            val = int(round(pct*total/100))
            return '{v:d}'.format(v=val)
        return my_format

    ax1.pie(ratios, autopct=autopct_format([pie1_1, pie1_2]), startangle=angle,
            labels=labels, colors=['#6d6e71', '#46bdf2'], explode=explode)

    # small pie chart parameters
    ratio_21 = pie2_1 / (pie2_1 + pie2_2)
    ratio_22 = pie2_2 / (pie2_1 + pie2_2)
    ratios = [ratio_21, ratio_22]
    labels = ['', '']
    # width = .2
    angle = -180 * ratios[1]
    radius = np.sqrt((pie2_1 + pie2_2) / (pie1_1 + pie1_2))
    ax2.pie(ratios, autopct=autopct_format([pie2_1, pie2_2]), startangle=angle,
            labels=labels, radius=radius, colors=['#84daff','#929497'],
            explode=explode, textprops={'size': 'smaller'})

    ax1.set_title('CORUM')
    ax2.set_title('BioPlex')
    # for now, no connection between pie charts
    # https://stackoverflow.com/a/72551175

    # small pie chart parameters
    ratio_31 = pie3_1 / (pie3_1 + pie3_2)
    ratio_32 = pie3_2 / (pie3_1 + pie3_2)
    ratios = [ratio_32, ratio_31]
    labels = ['', '']
    angle = -180 * ratios[0]
    radius = np.sqrt((pie3_1 + pie3_2) / (pie1_1 + pie1_2))
    ax3.pie(ratios, autopct=autopct_format([pie3_1, pie3_2]), startangle=angle,
            labels=labels, radius=radius, colors=['#d0d2d3', '#d3f3ff'],
            explode=explode, textprops={'size': 'smaller'})

    ax1.set_title('CORUM')
    ax2.set_title('BioPlex')
    ax3.set_title('GO')

    fig.tight_layout()
    plt.savefig(output_file, bbox_inches='tight')


if __name__ == '__main__':
    args = arg_parse()
    print('checking cytoscape connection')
    base_url = 'http://127.0.0.1:1234/v1'
    corum_file = 'allComplexes.txt'
    bioplex_edgelist_symbol_file = 'BioPlex_3_edge_list_gene_symbol.tsv'
    outfile = 'corum_complex_overlap_matrix.npy'
    outfile_filtered = 'filtered_corum_complexes_0.9.txt'
    ensembl_version = 'v104'
    release_date = 'may2021'
    out_dir = '.'
    id_mapping = None

    plot_pie_chart(args.enrich_results_file, os.path.join(out_dir, 'plot_enrichment_pie_chart.pdf'))

    print(p4c.cytoscape_ping(base_url=base_url))
    print(p4c.cytoscape_version_info(base_url=base_url))
    p4c.networks.delete_all_networks(base_url=base_url)
    p4c.styles.delete_all_visual_styles(base_url=base_url)

    g_name = 'plot_bioplex_corum_funmap'
    funmap_output_pdf = os.path.join(out_dir, f'{g_name}_network.pdf')
    funmap_output_cys = os.path.join(out_dir, f'{g_name}_network.cys')
    funmap_output_tsv = os.path.join(out_dir, f'{g_name}_network.tsv')

    if not os.path.exists(outfile_filtered):
        print(f'generating {outfile_filtered}')
        # id_mapping = load_id_map('hgnc_symbol', 'ensembl_gene_id',
        #                         ensembl_version=ensembl_version,
        #                         release_date=release_date)
        corum = load_corum_complexes(corum_file)
        corum['num_prot'] = corum.apply(lambda x: len(x['subunits']), axis=1)
        corum_filtered = corum.loc[corum['num_prot']>2,]
        # build a set of corum complexes that don't overlap more than 90% -> choose larger complex and discard
        # smaller one of both that overlap more than 90%
        # compute matrix with pairwise overlaps
        # delete smaller complex of largest overlaps and iterate until the largest overlap is below 90%
        if not os.path.exists(outfile):
            overlaps = np.zeros((corum_filtered.shape[0],corum_filtered.shape[0]),dtype=float)
            for i in range(corum_filtered.shape[0]):
                if i % 100 == 0:
                    print(i)
                for j in range(corum_filtered.shape[0]):
                    if i != j:
                        overlaps[i,j] = len(corum_filtered.iloc[i,2].intersection(corum_filtered.iloc[j,2]))/float(len(corum_filtered.iloc[i,2]))

            np.save(outfile,overlaps)
        else:
            print(f'loading overlaps from {outfile}')
            overlaps = np.load(outfile)

        cutoff = 0.9
        max_overlap = overlaps.max()
        ind = np.unravel_index(np.argmax(overlaps, axis=None), overlaps.shape)
        count = 0
        while max_overlap >= cutoff:
            overlaps[ind[0]] = -1
            max_overlap = overlaps.max()
            ind = np.unravel_index(np.argmax(overlaps, axis=None), overlaps.shape)
            count += 1
            if count % 100 == 0:
                print(count)

        target = open(outfile_filtered,'w')
        target.write('ComplexID\tComplexName\tComplexSubunits\n')
        for i in range(overlaps.shape[0]):
            if overlaps[i].min() >= 0:
                target.write(str(corum_filtered.iloc[i,0]) + '\t' + corum_filtered.iloc[i,1] + '\t' + '|'.join(list(corum_filtered.iloc[i,2])) + '\n')
        target.close()

    print(f'loading {outfile_filtered}')
    corum = pd.read_table(outfile_filtered)

    corum['num_pairs'] = [np.NaN for i in range(corum.shape[0])]
    corum['num_pairs'] = corum['num_pairs'].astype('Int64')
    corum['num_fm'] = [np.NaN for i in range(corum.shape[0])]
    corum['num_fm'] = corum['num_fm'].astype('Int64')
    corum['num_bp'] = [np.NaN for i in range(corum.shape[0])]
    corum['num_bp'] = corum['num_bp'].astype('Int64')
    corum['num_fm_bp'] = [np.NaN for i in range(corum.shape[0])]
    corum['num_fm_bp'] = corum['num_fm_bp'].astype('Int64')

    # if the file does not exist, create it
    if os.path.exists(bioplex_edgelist_symbol_file):
        print(f'loading {bioplex_edgelist_symbol_file}')
        bioplex_df = pd.read_table(bioplex_edgelist_symbol_file, header=None)
    else:
        print(f'generating {bioplex_edgelist_symbol_file}')
        bioplex_df = load_nw_bioplex(id_type='gene_symbol', out_file=bioplex_edgelist_symbol_file)

    bioplex = set(list(bioplex_df.itertuples(index=False,name=None)))

    # load funmap network
    funmap_df = pd.read_table(args.input_network_symbol, header=None)
    funmap_df.columns = ['gene_symbol_a', 'gene_symbol_b']
    funmap_df = funmap_df.set_index(funmap_df['gene_symbol_a'] + '_'
                                + funmap_df['gene_symbol_b'])
    funmap = set(list(funmap_df[['gene_symbol_a','gene_symbol_b']].itertuples(index=False,name=None)))

    nodes = []
    edges = []
    source_attr = []
    corum_ids = []
    symbols = []
    num_evidences = []

    for index,row in corum.iterrows():
        gene_symbols = row['ComplexSubunits'].split('|')
        complex_id = row['ComplexID']
        pairs = [tuple(sorted(list(t))) for t in itertools.combinations(row['ComplexSubunits'].split('|'),2)]
        corum.loc[corum['ComplexID'] == complex_id, 'num_pairs'] = len(pairs)
        pairs_set = set(pairs)
        funmap_pairs = pairs_set.intersection(funmap)
        corum.loc[corum['ComplexID'] == complex_id, 'num_fm'] = len(funmap_pairs)
        bioplex_pairs = pairs_set.intersection(bioplex)
        corum.loc[corum['ComplexID'] == complex_id, 'num_bp'] = len(bioplex_pairs)
        funmap_bioplex_pairs = funmap_pairs.intersection(bioplex_pairs)
        corum.loc[corum['ComplexID'] == complex_id, 'num_fm_bp'] = len(funmap_bioplex_pairs)

        if not (len(funmap_pairs) == 0 and len(bioplex_pairs) == 0):
            sources = ['corum' for i in range(len(pairs))]
            new_evidences = [0 for i in range(len(pairs))]
            for pair in funmap_pairs:
                sources[pairs.index(pair)] = 'funmap'
                new_evidences[pairs.index(pair)] = 0
            pairs_names = [(t[0] + '_' + str(complex_id), t[1] + '_' + str(complex_id)) for t in pairs]
            for pair in bioplex_pairs:
                if sources[pairs.index(pair)] == 'corum':
                    sources[pairs.index(pair)] = 'bioplex'
                else:
                    sources[pairs.index(pair)] = 'bioplex_funmap'

            symbols = symbols + gene_symbols
            corum_ids = corum_ids + [complex_id for i in range(len(symbols))]
            names = [gene_symbol + '_' + str(complex_id) for gene_symbol in gene_symbols]
            nodes = nodes + names
            edges = edges + pairs_names
            source_attr = source_attr + sources
            num_evidences = num_evidences + new_evidences

    corum.to_csv(funmap_output_tsv, sep='\t', index=False)
    # from IPython import embed; embed()
    print('creating cytoscape network...')

    g = igraph.Graph()
    g.add_vertices(nodes)
    g.add_edges(edges)
    g.es['origin'] = source_attr
    g.es['num_fm_evidences'] = num_evidences
    g.vs['complex_id'] = corum_ids
    g.vs['symbols'] = symbols

    network = p4c.networks.create_network_from_igraph(g, title=g_name, base_url=base_url)
    default_style = {
        'COMPOUND_NODE_PADDING': 10,
        'COMPOUND_NODE_SHAPE': 'ROUND_RECTANGLE',
        'NODE_BORDER_PAINT': '#000000',
        'NODE_BORDER_TRANSPARENCY': 255,
        'NODE_BORDER_WIDTH': 0.0,
        'NODE_FILL_COLOR': '#D4E5F4',
        'EDGE_BEND': '',
        'EDGE_WIDTH': 2.0,
        'EDGE_CURVED': True
    }
    node_label = p4c.style_mappings.map_visual_property('node label', 'symbols', 'p', base_url=base_url)
    edge_color = p4c.style_mappings.map_visual_property('edge stroke unselected paint',
                'origin', 'd', ['funmap', 'bioplex_funmap', 'bioplex', 'corum'],
                ['#FF0000', '#FF00FF', '#6699FF', '#FFFFFF'], base_url=base_url)
    edge_transparency = p4c.style_mappings.map_visual_property('edge transparency',
                'origin', 'd', ['funmap', 'bioplex_funmap', 'bioplex', 'corum'],
                [255, 255, 255, 0], base_url=base_url)
    style_name = 'corum_style'
    style_obj = p4c.styles.create_visual_style(style_name,
                defaults=default_style,
                mappings=[node_label, edge_color, edge_transparency], base_url=base_url)
    p4c.py4cytoscape_tuning.CATCHUP_NETWORK_TIMEOUT_SECS = 600
    p4c.layouts.layout_network(layout_name='force-directed', network=network, base_url=base_url)
    p4c.styles.set_visual_style(style_name, network=network, base_url=base_url)
    print('exporting image ...')
    p4c.network_views.fit_content(network=network, base_url=base_url)
    p4c.network_views.export_image(filename=funmap_output_pdf, type='PDF', overwrite_file=True, base_url=base_url)
    print('saving session ...')
    # saving seems to be working only with absolute path
    p4c.session.save_session(filename=os.path.abspath(funmap_output_cys), base_url=base_url)

    # from now on, we will use the full CORUM without filtering (with gene symbols)
    # load CORUM
    print('loading corum ...')
    corum_df = load_corum_complexes(corum_file)
    corum_edge_set = set()
    for index, row in corum_df.iterrows():
        gene_ids = list(row['subunits'])
        gene_ids = [gid.upper() for gid in gene_ids]
        pairs = set(tuple(sorted(list(t))) for t in itertools.combinations(gene_ids,2))
        corum_edge_set = corum_edge_set | pairs
    print('loading funmap ...')
    funmap_df = pd.read_csv(args.input_network_symbol, sep='\t', header=None)
    funmap_df.columns = ['gene_a', 'gene_b']
    funmap_edge_list = list(funmap_df[['gene_a','gene_b']].itertuples(index=False,name=None))
    funmap_edge_set = set(tuple(sorted([t[0].upper(), t[1].upper()])) for t in funmap_edge_list)
    print('loading bioplex ...')
    bioplex_df = pd.read_csv(bioplex_edgelist_symbol_file, sep='\t', header=None)
    bioplex_df.columns = ['gene_a', 'gene_b']
    bioplex_edge_list = list(bioplex_df[['gene_a','gene_b']].itertuples(index=False,name=None))
    bioplex_edge_set = set(tuple(sorted([t[0].upper(), t[1].upper()])) for t in bioplex_edge_list)

    sets = [funmap_edge_set, corum_edge_set, bioplex_edge_set]
    n_Abc = len(funmap_edge_set - corum_edge_set - bioplex_edge_set)
    n_aBc = len(corum_edge_set - funmap_edge_set - bioplex_edge_set)
    n_ABc = len(funmap_edge_set & corum_edge_set - bioplex_edge_set)
    n_abC = len(bioplex_edge_set - funmap_edge_set - corum_edge_set)
    n_AbC = len(funmap_edge_set & bioplex_edge_set -corum_edge_set)
    n_aBC = len(corum_edge_set & bioplex_edge_set - funmap_edge_set)
    n_ABC = len(funmap_edge_set & corum_edge_set & bioplex_edge_set)
    # for plotting, the numbers have different meaning
    print('n_Abc', 'n_aBc', 'n_ABc', 'n_abC', 'n_AbC', 'n_aBC', 'n_ABC')
    print(n_Abc, n_aBc, n_ABc, n_abC, n_AbC, n_aBC, n_ABC)
    print((n_ABc+n_ABC)/len(funmap_edge_set), (n_aBC+n_ABC)/len(bioplex_edge_set))
    print(len(funmap_edge_set), len(corum_edge_set), len(bioplex_edge_set))
    # plot venn diagram
    set_data = (n_Abc, n_aBc, n_ABc, n_abC, n_AbC, n_aBC, n_ABC)
    plot_venn(set_data, labels=('FunMap', 'CORUM', 'BioPlex'), filename='plot_funmap_corum_bioplex_venn.pdf')

    # plot comparions of funmap and bioplex ppis
    fig_data_file = 'plot_funmap_bioplex_dot_plot_data.tsv'
    corum_df_1 = corum_df.copy()
    corum_df_1['complex_size'] = corum_df_1.apply(lambda row: len(row['subunits']), axis=1)
    corum_df_1['num_pairs'] = corum_df_1.apply(lambda row: int(len(row['subunits'])*(len(row['subunits'])-1)/2), axis=1)
    num_fm_ppis_lst = []
    num_bp_ppis_lst = []
    num_fm_bp_ppis_lst = []

    for index, row in corum_df.iterrows():
        gene_ids = list(row['subunits'])
        gene_ids = [gid.upper() for gid in gene_ids]
        pairs = set(tuple(sorted(list(t))) for t in itertools.combinations(gene_ids,2))
        num_fm_ppis = len(funmap_edge_set & pairs)
        num_fm_ppis_lst.append(num_fm_ppis)
        num_bp_ppis = len(bioplex_edge_set & pairs)
        num_bp_ppis_lst.append(num_bp_ppis)
        num_fm_bp_ppis = len(funmap_edge_set & bioplex_edge_set & pairs)
        num_fm_bp_ppis_lst.append(num_fm_bp_ppis)

    corum_df_1['num_fm_ppis'] = num_fm_ppis_lst
    corum_df_1['num_bp_ppis'] = num_bp_ppis_lst
    corum_df_1['num_fm_bp_ppis'] = num_fm_bp_ppis_lst
    corum_df_1['subunits'] = corum_df_1.apply(lambda row: '|'.join(row['subunits']), axis=1)
    corum_df_1.to_csv(fig_data_file, sep='\t', index=False)

    funmap_pair_count = np.array(corum_df_1.loc[(corum_df_1['num_fm_ppis'] > 0) | (corum_df_1['num_bp_ppis'] > 0), 'num_fm_ppis']) + 1
    bp_pair_count = np.array(corum_df_1.loc[(corum_df_1['num_fm_ppis'] > 0) | (corum_df_1['num_bp_ppis'] > 0), 'num_bp_ppis']) + 1
    complex_size = np.array(corum_df_1.loc[(corum_df_1['num_fm_ppis'] > 0) | (corum_df_1['num_bp_ppis'] > 0), 'complex_size'])
    col_map = {'red': '#FF0000',
            'pink': '#FF00FF',
            'blue': '#6699FF'
    }
    fig_dotplot = plt.figure(figsize=(10,10))
    ax_bb = fig_dotplot.add_subplot(1, 1, 1)

    cur_size = len(funmap_pair_count[funmap_pair_count > bp_pair_count])
    scatter1 = ax_bb.scatter(bp_pair_count[funmap_pair_count > bp_pair_count],
                funmap_pair_count[funmap_pair_count > bp_pair_count],
                label=f'FunMap > BioPlex ({cur_size})',
                edgecolors='none', alpha=0.5, c=col_map['red'],
                s=complex_size[funmap_pair_count > bp_pair_count] * 10)
    cur_size = len(funmap_pair_count[funmap_pair_count == bp_pair_count])
    scatter2 = ax_bb.scatter(bp_pair_count[funmap_pair_count == bp_pair_count],
                funmap_pair_count[funmap_pair_count == bp_pair_count],
                label=f'FunMap = BioPlex ({cur_size})',
                edgecolors='none', alpha=0.5, c=col_map['pink'],
                 s=complex_size[funmap_pair_count == bp_pair_count] * 10)
    cur_size = len(funmap_pair_count[funmap_pair_count < bp_pair_count])
    scatter3 = ax_bb.scatter(bp_pair_count[funmap_pair_count < bp_pair_count],
                funmap_pair_count[funmap_pair_count < bp_pair_count],
                label=f'FunMap < BioPlex ({cur_size})',
                edgecolors='none', alpha=0.5, c=col_map['blue'],
                 s=complex_size[funmap_pair_count < bp_pair_count] * 10)

    lgnd1 = ax_bb.legend(loc='upper left', fontsize=24, frameon=True, handletextpad=0.2)
    for handle in lgnd1.legendHandles:
        handle.set_sizes([200.0])
    ax_bb.add_artist(lgnd1)
    max_val = max(max(funmap_pair_count), max(bp_pair_count))
    ax_bb.set_yscale('log')
    ax_bb.set_xscale('log')
    ax_bb.set_xlim(0.9, max_val+2000)
    ax_bb.set_ylim(0.9, max_val+2000)
    ax_bb.set_xlabel('BioPlex pairs count + 1')
    ax_bb.set_ylabel('Funmap pairs count + 1')
    ax_bb.set_aspect('equal')
    ax_bb.spines['top'].set_visible(False)
    ax_bb.spines['right'].set_visible(False)
    ax_bb.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_bb.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    handles1, labels1 = scatter1.legend_elements(prop='sizes', num=[10, 20, 50, 100], func=lambda x:x/10, alpha=1)
    lgnd2 = ax_bb.legend(handles1, labels1,
                        loc='lower right', labelspacing=0.7, fontsize=24,
                        title='Size', markerscale=1, frameon=True)
    for handle in lgnd2.legendHandles:
        handle.set_color('gray')

    ax_bb.add_artist(lgnd2)
    fig_dotplot.tight_layout()
    fig_dotplot.savefig(os.path.join('plot_funmap_bioplex_dot_plot.pdf'))

    fig_data_file = 'plot_funmap_corum_dot_plot_data.tsv'
    bioplex3_tsv = 'BioPlex3.tsv'
    bioplex_df = pd.read_csv(bioplex3_tsv, sep='\t')
    bioplex_df['gene_ids'] = [[] for _ in range(len(bioplex_df))]
    for index, row in bioplex_df.iterrows():
        id_set = row['genes'].split('|')
        if len(id_set) != 0:
            bioplex_df.loc[index, 'gene_ids'].extend(id_set)

    bioplex_df['complex_size'] = bioplex_df.apply(lambda row: len(row['gene_ids']), axis=1)

    num_fm_ppis_lst = []
    num_co_ppis_lst = []
    num_fm_co_ppis_lst = []
    num_pairs_lst = []

    for index, row in bioplex_df.iterrows():
        gene_ids = list(row['gene_ids'])
        gene_ids = [gid.upper() for gid in gene_ids]
        pairs = set(tuple(sorted(list(t)))
                    for t in itertools.combinations(gene_ids, 2))
        num_pairs = len(bioplex_edge_set & pairs)
        num_pairs_lst.append(num_pairs)
        real_pairs = bioplex_edge_set & pairs
        num_fm_ppis = len(funmap_edge_set & real_pairs)
        num_fm_ppis_lst.append(num_fm_ppis)
        num_co_ppis = len(corum_edge_set & real_pairs)
        num_co_ppis_lst.append(num_co_ppis)
        num_fm_co_ppis = len(funmap_edge_set & corum_edge_set & real_pairs)
        num_fm_co_ppis_lst.append(num_fm_co_ppis)

    bioplex_df['num_pairs'] = num_pairs_lst
    bioplex_df['num_fm_ppis'] = num_fm_ppis_lst
    bioplex_df['num_co_ppis'] = num_co_ppis_lst
    bioplex_df['num_fm_co_ppis'] = num_fm_co_ppis_lst
    bioplex_df['gene_ids'] = bioplex_df.apply(
        lambda row: '|'.join(row['gene_ids']), axis=1)
    bioplex_df.to_csv(fig_data_file, sep='\t', index=False)

    funmap_pair_count = np.array(bioplex_df.loc[(bioplex_df['num_fm_ppis'] > 0) | (
        bioplex_df['num_co_ppis'] > 0), 'num_fm_ppis']) + 1
    co_pair_count = np.array(bioplex_df.loc[(bioplex_df['num_fm_ppis'] > 0) | (
        bioplex_df['num_co_ppis'] > 0), 'num_co_ppis']) + 1
    complex_size = np.array(bioplex_df.loc[(bioplex_df['num_fm_ppis'] > 0) | (
        bioplex_df['num_co_ppis'] > 0), 'complex_size'])
    fig_b_1 = plt.figure(figsize=(10,10))
    ax_bb_1 = fig_b_1.add_subplot(1, 1, 1)
    cur_size = len(funmap_pair_count[funmap_pair_count > co_pair_count])
    scatter1 = ax_bb_1.scatter(co_pair_count[funmap_pair_count > co_pair_count],
                funmap_pair_count[funmap_pair_count > co_pair_count],
                label=f'FunMap > CORUM ({cur_size})',
                edgecolors='none', alpha=0.5, c=col_map['red'],
                s=complex_size[funmap_pair_count > co_pair_count] * 10)
    cur_size = len(funmap_pair_count[funmap_pair_count == co_pair_count])
    scatter2 = ax_bb_1.scatter(co_pair_count[funmap_pair_count == co_pair_count],
                funmap_pair_count[funmap_pair_count == co_pair_count],
                label=f'FunMap = CORUM ({cur_size})',
                edgecolors='none', alpha=0.5, c=col_map['pink'],
                 s=complex_size[funmap_pair_count == co_pair_count] * 10)
    cur_size = len(funmap_pair_count[funmap_pair_count < co_pair_count])
    scatter3 = ax_bb_1.scatter(co_pair_count[funmap_pair_count < co_pair_count],
                funmap_pair_count[funmap_pair_count < co_pair_count],
                label=f'FunMap < CORUM ({cur_size})',
                edgecolors='none', alpha=0.5, c=col_map['blue'],
                s=complex_size[funmap_pair_count < co_pair_count] * 10)

    lgnd1 = ax_bb_1.legend(loc='upper left', fontsize=24, frameon=True, handletextpad=0.2)
    for handle in lgnd1.legendHandles:
        handle.set_sizes([200.0])
    ax_bb_1.add_artist(lgnd1)
    max_val = max(max(funmap_pair_count), max(co_pair_count))
    ax_bb_1.set_yscale('log')
    ax_bb_1.set_xscale('log')
    ax_bb_1.set_xlim(0.75, max_val+2000)
    ax_bb_1.set_ylim(0.75, max_val+2000)
    ax_bb_1.set_xlabel('CORUM pairs count + 1')
    ax_bb_1.set_ylabel('Funmap pairs count + 1')
    ax_bb_1.set_aspect('equal')
    ax_bb_1.spines['top'].set_visible(False)
    ax_bb_1.spines['right'].set_visible(False)
    ax_bb_1.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax_bb_1.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())

    handles1, labels1 = scatter1.legend_elements(prop='sizes', num=[10, 20, 50, 100], func=lambda x:x/10, alpha=1)
    lgnd2 = ax_bb_1.legend(handles1, labels1,
                        loc='lower right', labelspacing=0.7, fontsize=24,
                        title='Complex size', markerscale=1, frameon=True)
    for handle in lgnd2.legendHandles:
        handle.set_color('gray')

    ax_bb_1.add_artist(lgnd2)
    fig_b_1.tight_layout()
    fig_b_1.savefig(os.path.join('plot_funmap_corum_dot_plot.pdf'))

    plot_ice_modules(args.input_network_symbol, args.ice_module_file,
                    base_url, default_style, 'plot_funmap_ice_modules.pdf')

    # note: if there are many cliques, this could take a while
    plot_ice_corum_overlap(args.enrich_results_file, default_style,
                        base_url, './ice_corum_overlap',
                        args.funmap_config_file,
                        args.data_file)
