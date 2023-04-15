#!/usr/bin/env python
import argparse
import sys
import os
from os import path as osp
import pandas as pd
import numpy as np
from scipy import stats
import csv
import json
from sklearn.metrics import make_scorer
from sklearn.model_selection import GridSearchCV
from xgboost import XGBRegressor
from sklearn.pipeline import Pipeline
import warnings
from sklearn.exceptions import ConvergenceWarning
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker as ticker
from statannotations.Annotator import Annotator
import seaborn as sns


warnings.simplefilter("ignore", category=ConvergenceWarning)

def arg_parse():
    parser = argparse.ArgumentParser(description='command line arguments.')
    parser.add_argument('-n', '--dataset-name', type=str, required=True,
                        help='name of the dataset')
    parser.add_argument('-m', '--module-info', type=str, required=True,
                        help='path to the module info file')
    parser.add_argument('-d', '--data-file', type=str, required=True,
                        help='path to the data file')
    parser.add_argument('-t', '--mutation-file', type=str, required=True,
                        help='path to the mutation file')
    return parser.parse_args()


def get_activity_scores(z_scores, module_dict):
    """ get the score for each module and clique for each sample
    """
    all_scores = None
    z_scores = z_scores.T
    for module in module_dict:
        cur_module = module_dict[module]
        common_genes = list(set(cur_module) & set(z_scores.columns))
        cur_score = z_scores.loc[:, common_genes].mean(axis=1).to_frame()
        cur_score.columns = [module]
        cur_module_score = cur_score

        if all_scores is None:
            all_scores = cur_module_score
        else:
            all_scores = pd.concat([all_scores, cur_module_score], axis=1)

    return all_scores


def prepare_features(mutation_file):
    df = pd.read_csv(mutation_file, sep='\t', index_col=0)
    df['Mutation_Frequency'] = df.sum(axis=1)
    df_sorted = df.sort_values('Mutation_Frequency', ascending=False)
    N = 100  # number of top mutated genes to use
    top_mutated_genes = df_sorted.head(N)
    top_mutated_genes_rows = df.loc[top_mutated_genes.index]
    top_mutated_genes_rows = top_mutated_genes_rows.drop(columns='Mutation_Frequency')

    return top_mutated_genes_rows.T


def corr_func(ground_truth, prediction):
    y_truth = np.reshape(ground_truth, -1)
    y_pred = np.reshape(prediction, -1)
    return stats.pearsonr(y_truth, y_pred)[0]


def xg_boost(features, target, n_jobs):
    idx = features.index.intersection(target.index)
    X = features.loc[idx, :]
    y = target.loc[idx, :]
    X_y = pd.concat([X, y], axis=1)
    # if more than half of samples contain missing values, skip
    percent_row_missing = X_y.isna().any(axis=1).sum() * 100 / X_y.shape[0]
    if percent_row_missing > 50.0:
        return None
    # remove samples without label, xgboost can take care of missing features
    X_y.dropna(axis=0, subset=[X_y.columns[-1]], inplace=True)
    X = X_y.iloc[:, :-1]
    y = X_y.iloc[:, -1]

    xgbr = XGBRegressor(n_jobs=n_jobs, random_state=1)
    pipe = Pipeline(steps=[('xgb', xgbr)])
    param_grid = {
        'xgb__learning_rate': [0.1, 0.2, 0.3, 0.4, 0.5],
        'xgb__n_estimators': [20, 50],
        'xgb__max_depth':[2, 3, 4]
    }

    my_scorer = make_scorer(corr_func, greater_is_better=True)
    search = GridSearchCV(pipe, param_grid, scoring=my_scorer, cv=5, n_jobs=1)
    search.fit(X, y)
    best_score = search.best_score_
    print(f'best score: {search.best_score_}')
    print(f'best params: {search.best_params_}')
    final_score = search.score(X, y)
    print(f'final score', final_score)
    sorted_features = None

    feature_importance = pd.Series(index=X.columns,
            data=search.best_estimator_.named_steps['xgb'].feature_importances_)
    sorted_features = feature_importance.sort_values(ascending=False)

    return {'best_score': best_score,
            'final_score': final_score,
            'feature_imp': sorted_features
            }


def predict_activity(features, targets, importance_file, best_score_file):
    """ train a model to predict module activity """
    # build one model for each target
    best_score = {'target':[], 'score':[]}
    feature_imp = {}
    for (idx, t) in enumerate(targets, start=1):
        print(f'{idx} of {targets.shape[1]}: {t}')
        sys.stdout.flush()
        results = xg_boost(features, targets[t].to_frame(), n_jobs=-1)
        if results is not None:
            best_score['target'].append(t)
            best_score['score'].append(results['best_score'])
            cur_feature_imp = results['feature_imp']
            if cur_feature_imp is not None:
                cur_feature_imp = cur_feature_imp.to_dict()
                for key in cur_feature_imp:
                    cur_feature_imp[key] = round(cur_feature_imp[key], 3)
                feature_imp[t] = cur_feature_imp

        with open(importance_file, 'w') as f:
            json.dump(feature_imp, f)

        best_score_df = pd.DataFrame(best_score)
        best_score_df.to_csv(best_score_file, index=False, sep='\t')


def plot_importance_barplot(output_dir, importance_file, top_n_bar):
    df = pd.read_csv(importance_file, sep='\t', index_col=0)
    df = df.sort_values('max', ascending=False)
    df = df.iloc[:top_n_bar, :]
    fig_name = os.path.splitext(importance_file)[0]
    fig_name = osp.join(output_dir, f'{fig_name}_barplot_{top_n_bar}.pdf')
    fig, ax = plt.subplots(figsize=(8, 4))
    ax.bar(df.index, df['max'], color='#E72F52', label='max')
    ax.bar(df.index, -df['total'], color='#0D95D0', label='total')
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 8))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.yaxis.grid(color = 'gainsboro', linewidth=2, linestyle = 'dotted')
    ax.xaxis.set_ticks(range(len(df.index)))
    ax.set_xticklabels(list(df.index), rotation=90)
    ax.set_xlim(-1, len(df.index))
    ax.set_axisbelow(True)
    ax.set_ylabel('Importance')
    def positive(x, pos):
        return f'{abs(x)}'

    formatter = FuncFormatter(positive)
    ax.yaxis.set_major_formatter(formatter)
    plt.legend()
    plt.tight_layout()
    fig.savefig(fig_name,  bbox_inches = 'tight', pad_inches = 0.1)


def plot_mutation_activity(out_dir, dataset_name, importance_file, feature_file,
                            activity_score_file):
    feature_df= pd.read_csv(feature_file, sep='\t', index_col=0)
    activity_scores = pd.read_csv(activity_score_file, sep='\t', index_col=0)

    imp_df = pd.read_csv(importance_file, sep='\t', index_col=0)
    # remove the last 2 columns ('max', 'total')
    imp_df = imp_df.iloc[:, :-2]
    # find the row and column with the top 10 largest values
    top_values = imp_df.stack().nlargest(10)
    row_col = top_values.index.tolist()
    gene_module = [(rc[0], rc[1]) for rc in row_col]
    for i in gene_module:
        print(imp_df.loc[i[0], i[1]])

    samples_dict = {}
    ct = dataset_name
    for idx, (gene, module) in enumerate(gene_module, start=1):
        all_data = {}
        cur_data = feature_df.loc[:, gene]
        mut = cur_data[cur_data == 1]
        no_mut = cur_data[cur_data == 0]
        # require both (mutated , not-mutated) at least 3 samples
        if len(mut) < 3 or len(no_mut) < 3:
            print(f'not enough samples for {gene} {module} (mutated: {len(mut)}, not-mutated: {len(no_mut)})')
            continue
        all_data[ct] = (len(no_mut), len(mut))
        samples_dict[ct] = feature_df.index.tolist()

        h = 5
        w = 3
        fig, ax = plt.subplots(1, 1, figsize=(w, h))
        fig_name = os.path.join(out_dir, f'mutation_activity_{idx}_{gene}__{module}.pdf')
        print(fig_name)
        for (i, ct) in enumerate(all_data):
            plot_data = pd.DataFrame({'mutation':feature_df.loc[samples_dict[ct], gene],
                                    'score':activity_scores.loc[samples_dict[ct], module]})
            map_dict = {0: 'wt', 1:'mut'}
            plot_data['mutation'] = plot_data['mutation'].map(map_dict)
            plotting = {
                'data': plot_data,
                'x': 'mutation',
                'y': 'score',
                'order': ['wt', 'mut']
            }
            my_pal = {'wt': '#a1d76a', 'mut': '#e9a3c9'}
            sns.boxplot(ax=ax, width=0.5, palette=my_pal, linewidth=0.8,
                        **plotting)
            ax.set_title(ct)
            ax.set_xticklabels([f'wt\n({all_data[ct][0]})',
                                f'mut\n({all_data[ct][1]})'])
            ax.set_xlabel(gene)
            ax.yaxis.grid(color='gainsboro', linestyle='dotted')
            ax.set_axisbelow(True)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_visible(False)
            ax.spines['bottom'].set_position(('outward', 5))
            if i == 0:
                ax.set_ylabel(f'{module} activity score')
            else:
                ax.set_ylabel('')
                ax.set_yticklabels([])
                ax.tick_params(axis='y', top=False, bottom=False, left=False, right=False,
                    labelleft=False, labelbottom=False)
            annotator = Annotator(ax=ax, pairs=[('wt', 'mut')], data=plot_data,
                                x='mutation', y='score')
            annotator.configure(test='t-test_ind', text_format='star')
            annotator.apply_and_annotate()
        plt.tight_layout()
        fig.savefig(fig_name,  bbox_inches = 'tight', pad_inches = 0.1)
        plt.close(fig)


def plot_importance_dotplot(out_dir, importance_file, score_file, threshold_cc, top_n=20):
    scores = pd.read_csv(score_file, sep='\t', index_col=0)
    score_sorted = scores.sort_values(by='score', ascending=False)
    score_selected = score_sorted[score_sorted['score'] > threshold_cc]
    if len(score_selected) > top_n:
        score_selected = score_sorted.iloc[:top_n, :]
    sorted_modules = score_selected.index.tolist()
    df = pd.read_csv(importance_file, sep='\t', index_col=0)
    df = df.sort_values('max', ascending=False)
    df = df.iloc[:top_n, :-2]
    max_val = round(df.max().max(), 1)
    df = df[sorted_modules]
    fig_name = os.path.splitext(importance_file)[0]
    fig_name = os.path.join(out_dir, f'{fig_name}_dotplot_{top_n}.pdf')
    n_row = top_n
    n_col = df.shape[1]
    width = 8
    height_mul = 0.6
    height = np.ceil(n_row * width / n_col) * height_mul
    fig, axs = plt.subplots(nrows=2, ncols=2,
                        sharex=True, gridspec_kw={'height_ratios': [2, 0.5], 'width_ratios': [1, 0.1]})
    ax, ax_3, ax_2, ax_4 = axs.ravel()
    fig.set_size_inches(width*1.01, 1.25*height)
    fig.subplots_adjust(hspace=0.5)
    ax.set_axisbelow(True)
    ax.xaxis.grid(color = 'gainsboro', linestyle = 'dotted')
    ax.yaxis.grid(color = 'gainsboro', linestyle = 'dotted')
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    edge_colors = None
    size_mul = 500
    for i in range(n_row):
        x = np.array(range(n_col))
        y = np.array([i] * n_col)
        v = np.array(df.iloc[n_row-1-i, :].to_list())
        sz = np.abs(np.array(df.iloc[n_row-1-i, :].to_list()))
        sc = ax.scatter(x, y, c=v, edgecolors=edge_colors, vmin=0,
                        vmax=max_val, cmap='Reds', s=sz*size_mul)

    ax.xaxis.set_ticks(np.arange(n_col))
    ax.yaxis.set_ticks(np.arange(n_row))
    xtick_lab = df.columns.to_list()
    ax.set_xticklabels(xtick_lab, rotation=90)
    features = df.index.to_list()
    features.reverse()
    ax.set_yticklabels(features)
    cb = plt.colorbar(sc, shrink=0.2, ax=ax_3)
    cb.outline.set_edgecolor('white')
    cb.ax.set_ylim(0, max_val)
    cb_locator = ticker.MultipleLocator(base=0.2)  # Set tick locations every 0.2
    cb.ax.xaxis.set_major_locator(cb_locator)
    L = ax_3.legend(*sc.legend_elements('sizes', num=4), title='feature\nimportance', borderpad=1,
                    labelspacing=1, loc='upper left', frameon=False)
    text_len = len(L.get_texts())
    real_numbers = []
    for text in L.get_texts():
        cur_label = text.get_text()  # Get the text from the Text object
        number = float(cur_label.strip('$\\mathdefault{}$'))  # Extract number by removing leading and trailing characters
        real_numbers.append(round(number/size_mul, 3))
    for i in range(text_len):
        L.get_texts()[i].set_text(real_numbers[i])

    leg = ax_3.get_legend()
    for i in range(len(leg.legendHandles)):
        leg.legendHandles[i].set_color('gray')

    ax_2.spines['bottom'].set_visible(False)
    ax_2.spines['right'].set_visible(False)
    ax_2.set_axisbelow(True)
    bar_color = '#0D95D0'
    ax_2.yaxis.grid(color = 'gainsboro', linestyle = 'dotted')
    ax_2.bar(score_selected.index, score_selected['score'], color=bar_color)
    ax_2.set_ylim(1, 0)
    ax_2.set_ylabel('Prediction\nPerformance (PCC)')
    ax_2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    ax.tick_params(axis='x', which='both', labelbottom=True)
    ax_2.tick_params(axis='x', which='both', top=False, bottom=False, labelbottom=False)
    ax_3.axis('off')
    ax_4.axis('off')

    plt.tight_layout()
    fig.savefig(fig_name,  bbox_inches = 'tight', pad_inches = 0.1)
    return sorted_modules


def plot_performance_barplot(out_dir, score_file):
    bar_color = '#0D95D0'
    scores = pd.read_csv(score_file, sep='\t', index_col=0)
    df_sorted = scores.sort_values(by='score', ascending=False)
    heights_sorted = df_sorted['score'].reset_index(drop=True)
    num_bars = len(heights_sorted)
    fig, ax = plt.subplots(figsize=(12,6))
    # Create bar plot
    ax.bar(np.arange(num_bars), heights_sorted, color=bar_color)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 8))
    ax.spines['bottom'].set_position(('outward', 5))
    ax.yaxis.grid(color = 'gainsboro', linewidth=2, linestyle = 'dotted')
    ax.set_axisbelow(True)
    ax.set_ylabel('Performance score (PCC)')
    ax.set_xlabel('Module')
    ax.xaxis.set_ticks(np.arange(num_bars))
    ax.xaxis.set_ticklabels([])
    plt.tight_layout()
    fig.savefig(os.path.join(out_dir, 'performance_barplot.pdf'),
                bbox_inches = 'tight', pad_inches = 0.1)


if __name__ == '__main__':
    args = arg_parse()
    output_dir = '.'
    module_info = args.module_info
    # mutation file: rows are genes, columns are samples, values are 0/1
    mutation_file = args.mutation_file
    data_file = args.data_file
    # only select cc > 0.25 for plotting
    threshold_cc = 0.25

    print('preparing input data ...')
    data_df = pd.read_csv(data_file, delimiter='\t', na_values=['NA', ''], index_col=0)
    z_scores = stats.zscore(data_df, axis=1, nan_policy='omit')
    z_scores = pd.DataFrame(z_scores, index=data_df.index, columns=data_df.columns)

    module_dict = {}
    with open(module_info, 'r', newline='') as file:
        reader = csv.reader(file, delimiter='\t')
        for row in reader:
            key = row[0]
            values = row[1:]
            module_dict[key] = values

    # activity_scores: rows are samples, columns are modules, values are
    # activity scores for each module in each sample
    activity_scores = get_activity_scores(z_scores, module_dict)
    activity_score_file = osp.join(output_dir, 'activity_scores.tsv')
    activity_scores.to_csv(activity_score_file, sep='\t', index=True)
    # features: rows are samples, columns are genes, values are 0/1
    features = prepare_features(mutation_file)
    feature_file = osp.join(output_dir, 'features.tsv')
    features.to_csv(feature_file, sep='\t', index=True)
    print('data preparation done.')

    importance_file = osp.join(output_dir, 'feature_importance.json')
    best_score_file = osp.join(output_dir, 'best_score.tsv')
    print('predict activity: model training ...')
    predict_activity(features, activity_scores, importance_file, best_score_file)
    print('model training completed.')
    best_score = pd.read_csv(best_score_file, sep='\t')
    with open(importance_file, 'r') as f:
        feature_imp = json.load(f)
    above_threshold = best_score[best_score['score'] > threshold_cc]

    if above_threshold.shape[0] > 20:
        best_score = best_score.sort_values(by=['score'], ascending=False)
        selected = best_score.iloc[:20, :]
    else:
        selected = above_threshold

    selected_modules = selected['target'].to_list()
    combined = None
    for i in feature_imp:
        if i not in selected_modules:
            continue
        cur_df = pd.DataFrame.from_dict(feature_imp[i], orient='index', columns=[i])
        if combined is None:
            combined = cur_df
        else:
            combined = pd.merge(combined, cur_df, how='outer',
                                right_index=True, left_index=True)
    combined.fillna(0, inplace=True)
    abs_total = combined.abs().sum(axis=1)
    abs_max = combined.abs().max(axis=1)
    combined['total'] = abs_total
    combined['max'] = abs_max
    importance_file_selected = f'feature_importance_selected.tsv'
    combined.to_csv(osp.join(output_dir, importance_file_selected), sep='\t')

    plot_performance_barplot(output_dir, best_score_file)
    plot_importance_barplot(output_dir, importance_file_selected, top_n_bar=50)
    selected_modules = plot_importance_dotplot(output_dir, importance_file_selected, best_score_file,
                            threshold_cc, top_n=20)
    with open(osp.join(output_dir, 'selected_modules.txt'), 'w') as f:
        for item in selected_modules:
            f.write(item + '\n')
    plot_mutation_activity(output_dir, args.dataset_name, importance_file_selected, feature_file,
                        activity_score_file)
