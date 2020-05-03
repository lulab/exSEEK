from collections import defaultdict
import pandas as pd
import h5py
import numpy as np
from copy import deepcopy
from tqdm import tqdm
import os

has_feature_stability = True
if snakemake.wildcards.cross_validation == 'cross_validation_diffexp':
    def parse_dirname(dirname):
        c = dirname.split('/')
        data = {}
        data['classifier'], data['n_features'], data['diffexp_method'], data['fold_change_direction'] = c[-1].split('.')
        data['compare_group'] = c[-2]
        data['count_method'] = c[-3]
        return data
elif snakemake.wildcards.cross_validation == 'cross_validation':
    def parse_dirname(dirname):
        c = dirname.split('/')
        data = {}
        data['classifier'], data['n_features'], data['selector'], data['fold_change_direction'] = c[-1].split('.')
        data['compare_group'] = c[-2]
        data['filter_method'], data['imputation'], data['normalization'], data['batch_removal'], data['count_method'] = c[-3].split('.')
        data['preprocess_method'] = '.'.join((data['filter_method'], data['imputation'], data['normalization'], data['batch_removal']))
        return data
elif snakemake.wildcards.cross_validation == 'evaluate_features':
    has_feature_stability = False
    def parse_dirname(dirname):
        c = dirname.split('/')
        data = {}
        data['classifier'] = c[-1]
        data['filter_method'], data['imputation'], data['normalization'], data['batch_removal'], data['count_method'] = c[-2].split('.')
        data['preprocess_method'] = '.'.join((data['filter_method'], data['imputation'], data['normalization'], data['batch_removal']))
        data['feature_set'] = c[-3]
        data['compare_group'] = c[-4]
        return data
else:
    raise ValueError('unknown cross_validation directory: {}'.format(snakemake.wildcards.cross_validation))

def feature_stability(X):
    '''Feature stability based on Kuncheva index
    
    Parameters:
    ----------
    
    X: array-like, shape (n_resample, n_features)
        Boolean matrix indicating features selected in each resamping run
    
    Returns:
    --------
    
    stability: float
        Feature stability
    '''
    # signature size
    s = X.sum(axis=1).mean()
    # total number of features
    N = X.shape[1]
    # number of resamping runs
    k = X.shape[0]
    # number of common signatures
    X = X.astype(np.int32)
    r = X.dot(X.T)[np.tril_indices(k, k=-1)]
    # Kuncheva index
    KI = (r - ((s**2)/N))/(s - ((s**2)/N))
    # pairwise stability score
    stability = 2*KI.sum()/(k*(k - 1))
    return stability

# read selected methods for each clustering score
# clustering_score[count_method][preprocess_method] = score_name
'''
clustering_score_names = defaultdict(dict)
for filename in snakemake.input.selected_methods:
    c = filename.split('/')
    count_method = c[-2]
    score_name = c[-3]
    with open(filename, 'r') as f:
        preprocess_method = f.readline().strip()
    if preprocess_method not in clustering_score_names[count_method]:
        clustering_score_names[count_method][preprocess_method] = []
    clustering_score_names[count_method][preprocess_method].append(score_name)
'''

columns = defaultdict(list)
summary = defaultdict(list)
for input_dir in tqdm(snakemake.input.input_dir, unit='directory'):
    metadata = parse_dirname(input_dir)
    # add clustering score that select the preprocess_method
    #if 'preprocess_method' in metadata:
        #names = clustering_score_names[metadata['count_method']][metadata['preprocess_method']]
    #else:
    #    names = ['null']
    #for clustering_score_name in names:
    #metadata = deepcopy(metadata)
    #metadata['clustering_score_name'] = clustering_score_name
    # feature_stability
    if has_feature_stability:
        with h5py.File(os.path.join(input_dir, 'cross_validation.h5'), 'r') as f:
            feature_selection_matrix = f['feature_selection'][:]
        record = deepcopy(metadata)
        record['feature_stability'] = feature_stability(feature_selection_matrix)
        summary['feature_stability'].append(record)
        if not columns['feature_stability']:
            columns['feature_stability'] = list(record.keys())
    # metrics
    for subset in ('train', 'test'):
        df_metrics = pd.read_table(os.path.join(input_dir, 'metrics.%s.txt'%subset), sep='\t')
        df_metadata = pd.DataFrame(index=np.arange(df_metrics.shape[0]))
        for key, val in metadata.items():
            df_metadata[key] = val
        columns['metrics.%s'%subset] = list(metadata.keys()) + df_metrics.columns.tolist()
        summary['metrics.%s'%subset].append(pd.concat([df_metadata, df_metrics], axis=1))
        
summary['metrics.train'] = pd.concat(summary['metrics.train'], axis=0)
summary['metrics.test'] = pd.concat(summary['metrics.test'], axis=0)
summary['metrics.train'] = summary['metrics.train'].reindex(columns=columns['metrics.train'])
summary['metrics.test'] = summary['metrics.test'].reindex(columns=columns['metrics.test'])
summary['metrics.train'].to_csv(snakemake.output.metrics_train, sep='\t', header=True, index=False, na_rep='NA')
summary['metrics.test'].to_csv(snakemake.output.metrics_test, sep='\t', header=True, index=False, na_rep='NA')
if has_feature_stability:
    summary['feature_stability'] = pd.DataFrame.from_records(summary['feature_stability'])
    summary['feature_stability'] = summary['feature_stability'].reindex(columns=columns['feature_stability'])
    summary['feature_stability'].to_csv(snakemake.output.feature_stability, sep='\t', header=True, index=False, na_rep='NA')