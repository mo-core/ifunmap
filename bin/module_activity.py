#!/usr/bin/env python
import os
import csv
import argparse
import pandas as pd
import numpy as np
from io import StringIO
import scipy.stats as stats
from scipy.stats import zscore, mannwhitneyu, f_oneway
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from lifelines import CoxPHFitter

plt.rcParams['xtick.labelsize'] = 6
plt.rcParams['ytick.labelsize'] = 6

def get_module_list_from_nsm(nsm_file):
    """
    get the module list from the nsm file, each module is a list of genes symbols
    nsm_file: the nsm file from NetSAM
    """
    lines = []
    with open(nsm_file,'r') as f:
        flag = False
        for line in f:
            if line.strip().startswith('##  Ruler file'):
                flag = True
                continue
            elif line.strip().startswith('##'):
                flag = False
            if flag:
                lines.append(line)

    ruler_df = pd.read_csv(StringIO('\n'.join(lines)), sep='\t')
    lines = []
    with open(nsm_file,'r') as f:
        flag = False
        for line in f:
            if line.strip().startswith('##  HMI file'):
                flag = True
                continue
            elif line.strip().startswith('##'):
                flag = False
            if flag:
                lines.append(line)

    hmi_df = pd.read_csv(StringIO('\n'.join(lines)), sep='\t')

    module_list = []
    module_start_end = []

    for _, row in hmi_df.iterrows():
        if row['level'] == 0:
            continue
        else:
            cur_module = {}
            cur_module['level'] = row['level']
            cur_module['order'] = row['order']
            cur_module['name'] = row['name']
            cur_module['name'] = cur_module['name'].replace('Level_', 'L')
            cur_module['name'] = cur_module['name'].replace('Module_', 'M')
            cur_module['start'] = row['start']
            cur_module['end'] = row['end']
            start = row['start'] - 1
            end = row['end'] - 1
            cur_module['genes'] = ruler_df.loc[start:end, 'node_name'].to_list()
            module_list.append(cur_module)
            cur_module_start_end = {}
            cur_module_start_end['level'] = row['level']
            cur_module_start_end['order'] = row['order']
            cur_module_start_end['name'] = cur_module['name']
            cur_module_start_end['start'] = row['start']
            cur_module_start_end['end'] = row['end']
            cur_module_start_end['size'] = row['end'] - row['start'] + 1
            module_start_end.append(cur_module_start_end)

    return (module_list, module_start_end)


def plot_activity_score(results_file, output_file):
    data = pd.read_csv(results_file, delimiter='\t')
    y_data = data.iloc[:, 1:]

    abs_data = np.abs(data.iloc[:, 1:])
    abs_data = np.nan_to_num(abs_data)
    top_5 = np.unravel_index(np.argsort(abs_data.ravel())[-5:], abs_data.shape)

    # create figure and axes
    fig, (ax_pos, ax_neg) = plt.subplots(nrows=2, ncols=1,
                                        sharex=True, gridspec_kw={'height_ratios': [3, 3]})
    fig.subplots_adjust(hspace=0.5)
    max_val = y_data.max().max()
    min_val = y_data.min().min()
    # find the absolute max value
    max_abs = max(abs(min_val), abs(max_val))
    y_lim = max_abs + 0.5
    # round y_lim to the nearest integer
    y_lim = int(round(y_lim))

    # plot positive values in the first subplot
    for i, col in enumerate(y_data.columns):
        cur_data = y_data[col]
        cur_data_pos = cur_data[cur_data > 0]
        age_jitter_pos = np.random.normal(scale=0.05, size=len(cur_data_pos))
        ax_pos.plot(age_jitter_pos+i, cur_data_pos, '.',  markeredgecolor='none', alpha=0.5, label=col)
        cur_data_neg = cur_data[cur_data < 0]
        age_jitter_neg = np.random.normal(scale=0.05, size=len(cur_data_neg))
        ax_neg.plot(age_jitter_neg+i, cur_data_neg, '.', alpha=0.5,  markeredgecolor='none', label=col)

    ax_pos.set_ylim(0, y_lim)
    # ax_pos.set_yticks(np.arange(0, y_lim+1, 2))
    major_yticks = np.arange(0, y_lim+1, 2)
    ax_pos.yaxis.set_major_locator(ticker.MultipleLocator(2))
    ax_pos.set_yticks(major_yticks)
    # ax_pos.set_yticks(ticker.MaxNLocator(nbins=10).tick_values(0, y_lim))
    ax_pos.tick_params(axis='x', labelbottom=True)
    ax_pos.axhline(y=2, color='gray', linestyle='--', linewidth=0.5)

    ax_neg.set_ylim(-y_lim, 0)
    # ax_neg.set_yticks(ticker.MaxNLocator(nbins=10).tick_values(-y_lim, 0))
    ax_neg.set_yticks(np.arange(-y_lim, 1, 2))
    ax_neg.yaxis.set_major_locator(ticker.MultipleLocator(2))
    # ax_neg.set_yticks(np.arange(-y_lim, 0, 1))
    # ax_neg.set_ylabel('Negative values')
    # hide x-axis tick labels for the bottom subplot
    ax_neg.tick_params(axis='x', labelbottom=False)
    ax_neg.axhline(y=-2, color='gray', linestyle='--', linewidth=0.5)
    ax_neg.set_xticklabels([])  # set the xticklabels to an empty list

    # add x-axis tick labels
    ax_pos.set_xticks(range(len(y_data.columns)))
    ax_pos.set_xticklabels(y_data.columns, rotation=90)

    ax_pos.tick_params(axis='x', which='both', labelbottom=True)
    ax_neg.tick_params(axis='x', which='both', top=False, bottom=False, labelbottom=False)

    ax_pos.spines['top'].set_visible(False)
    ax_pos.spines['right'].set_visible(False)
    ax_neg.spines['bottom'].set_visible(False)
    ax_neg.spines['right'].set_visible(False)
    fig.text(0, 0.5, 'signed -log10(p-value)', ha='center', va='center', rotation='vertical')

    jitter_amount = 0.2
    modules = data.iloc[:, 0]
    for i in range(len(top_5[0])):
        row, col = top_5[0][i], top_5[1][i]
        value = data.iloc[row, col + 1]
        module_name = modules[row]
        annotation = f'{module_name}\n{value:.2f}'
        jitter_x = np.random.uniform(-jitter_amount, jitter_amount)
        jitter_y = np.random.uniform(-jitter_amount, jitter_amount)
        if value > 0:
            ax_pos.annotate(annotation, xy=(col, value), xytext=(col+jitter_x, value+jitter_y),
                            arrowprops=dict(facecolor='red', edgecolor='red', arrowstyle='->'),
                            color='black', ha='center', va='center')
        else:
            ax_neg.annotate(annotation, xy=(col, value), xytext=(col+jitter_x, value+jitter_y),
                            arrowprops=dict(facecolor='red', edgecolor='red', arrowstyle='->'),
                            color='black', ha='center', va='center')
    fig.savefig(output_file, bbox_inches='tight')



def arg_parse():
    parser = argparse.ArgumentParser(description='command line arguments.')
    parser.add_argument('-i', '--ice-clique', type=str, required=True,
                        help='path to the clique file from ICE.')
    parser.add_argument('-n', '--netsam-nsm', type=str, required=True,
                        help='path to the nsm file from NetSAM.')
    parser.add_argument('-d', '--data-file', type=str, required=True,
                        help='path to the data matrix file')
    parser.add_argument('-t', '--tsi-file', type=str, required=True,
                        help='path to the tsi file')
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parse()
    module_dict = {}
    with open(args.ice_clique, 'r') as fh:
        for line in fh:
            line = line.strip()
            if line.startswith('#'):
                continue
            genes = line.split('\t')
            # the first element is the module name, the rest are genes
            module_name = genes[0]
            module_dict[module_name] = genes[1:]

    module_list, _ = get_module_list_from_nsm(args.netsam_nsm)
    for (i, module) in enumerate(module_list):
            module_dict[module['name']] = module['genes']

    # save the module dict to npy file
    np.save('module_info.npy', module_dict)
    with open('module_info.tsv', 'w', newline='') as file:
        writer = csv.writer(file, delimiter='\t')
        for key, values in module_dict.items():
            row = [key] + values
            writer.writerow(row)

    # read the tsi file
    tsi_df = pd.read_csv(args.tsi_file, delimiter='\t', skiprows=[1],
                        na_values=['NA', ''], index_col=0)
    # metadata includes the data type of each column
    with open(args.tsi_file) as f:
        metadata_row = f.readlines()[1].strip().split('\t')[1:]
    metadata = dict(zip(tsi_df.columns, metadata_row))

    # read the data file
    data_df = pd.read_csv(args.data_file, delimiter='\t', na_values=['NA', ''], index_col=0)
    z_scores = zscore(data_df, axis=1, nan_policy='omit')
    z_scores = pd.DataFrame(z_scores, index=data_df.index, columns=data_df.columns)

    all_results = None
    for feature in metadata:
        print(feature)
        results_df = pd.DataFrame(columns=['module', feature])
        for module in module_dict:
            genes = module_dict[module]
            # filter out genes that are not in the data file
            genes = [gene for gene in genes if gene in data_df.index]
            if len(genes) == 0:
                results_df = results_df.append({'module': module, feature: np.nan}, ignore_index=True)
                continue
            if metadata[feature] == 'CON':
                module_zscore = z_scores.loc[genes,:]
                module_activity = module_zscore.mean(axis='index').to_frame(name='activity')
                feature_df = tsi_df.loc[:, [feature]]
                # merge based on the index
                merged = pd.merge(module_activity, feature_df, left_index=True, right_index=True)
                # remove rows with na
                merged.dropna(inplace=True)
                cc_val, p_value = stats.spearmanr(merged.iloc[:,0], merged.iloc[:,1])
                score = np.sign(cc_val) * (-np.log10(p_value))
                results_df = results_df.append({'module': module, feature: score}, ignore_index=True)
            elif metadata[feature] == 'BIN':
                feature_values = tsi_df.loc[:, feature].unique()
                # if there is na in it then remove it
                feature_values = [value for value in feature_values if not pd.isna(value)]
                feature_values.sort()
                if len(feature_values) != 2:
                    raise ValueError(f'feature {feature} is not binary')
                # find out sample names with feature value 0 and 1
                sample_1 = tsi_df[tsi_df[feature] == feature_values[0]]
                sample_1 = sample_1.index
                sample_2 = tsi_df[tsi_df[feature] == feature_values[1]]
                sample_2 = sample_2.index
                # get the z-scores of the samples
                data_1 = z_scores.loc[genes, sample_1.to_list()].mean(axis='index').values
                data_2 = z_scores.loc[genes, sample_2.to_list()].mean(axis='index').values
                statistic, p_value = mannwhitneyu(data_1, data_2)
                direction = np.sign(np.mean(data_2) - np.mean(data_1))
                signed_log_p_value = direction * (-np.log10(p_value))
                results_df = results_df.append({'module': module, feature: signed_log_p_value},
                                            ignore_index=True)
            elif metadata[feature] == 'CAT':
                # we will use one-way ANOVA to test the significance of the module activity
                num_categories = len(tsi_df.loc[:, feature].unique())
                num_categories = num_categories - 1 if tsi_df.loc[:, feature].isna().sum() > 0 else num_categories
                if num_categories < 2:
                    raise ValueError(f'feature {feature} does not have enough categories')
                if num_categories == 2: # should be BIN
                    raise ValueError(f'feature {feature} is binary')
                module_zscore = z_scores.loc[genes,:]
                module_activity = module_zscore.mean(axis='index').to_frame(name='activity')
                feature_df = tsi_df.loc[:, [feature]]
                merged = pd.merge(module_activity, feature_df, left_index=True, right_index=True)
                merged = merged.dropna()
                categories = merged.loc[:, feature].unique()
                groups = [merged[merged[feature] == category]['activity'] for category in categories]
                f_statistic, p_value = f_oneway(*groups)
                score = np.sign(f_statistic) * (-np.log10(p_value))
                results_df = results_df.append({'module': module, feature: score}, ignore_index=True)
            elif metadata[feature] == 'SUR':
                # we will perform survival analysis
                module_zscore = z_scores.loc[genes,:]
                module_activity = module_zscore.mean(axis='index').to_frame(name='activity')
                feature_df = tsi_df.loc[:, [feature]]
                merged = pd.merge(module_activity, feature_df, left_index=True, right_index=True)
                # the survival data is in the format of "time,event", it is
                # possible that the time or event is missing
                # build the model
                merged[['time_to_event', 'event_occurred']] = merged[feature].str.split(',', expand=True)
                merged['time_to_event'] = pd.to_numeric(merged['time_to_event'], errors='coerce')
                merged['event_occurred'] = pd.to_numeric(merged['event_occurred'], errors='coerce')
                merged.drop(columns=[feature], inplace=True)
                merged.dropna(inplace=True)
                cph = CoxPHFitter()
                # fit the model with the continuous predictor variable
                cph.fit(merged, duration_col='time_to_event', event_col='event_occurred',
                        formula='activity')
                # get the p-value of the predictor variable
                p_val = cph.summary.loc['activity', 'p']
                # get the stats of the predictor variable
                stat_val = cph.summary.loc['activity', 'coef']
                score = np.sign(stat_val) * (-np.log10(p_val))
                results_df = results_df.append({'module': module, feature: score}, ignore_index=True)

        if all_results is None:
            all_results = results_df
        else:
            all_results = pd.merge(all_results, results_df, on='module')

    results_file = 'module_activity_scores.tsv'
    fig_name = os.path.splitext(results_file)[0]
    fig_name = f'{fig_name}.pdf'
    all_results.to_csv(results_file, sep='\t', index=False)
    plot_activity_score(results_file, fig_name)
