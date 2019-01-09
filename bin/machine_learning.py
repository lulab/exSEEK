#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

import json
from abc import ABC, abstractmethod
from tqdm import tqdm

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f


def parse_params(s):
    if s:
        params = s
        try:
            params = json.loads(s)
        except json.JSONDecodeError as e:
            pass
        finally:
            return params
    else:
        return {}

def make_estimator(
    feature_names=None,
    zero_fraction_filter=None,
    fold_change_filter=None,
    rpkm_filter=None,
    log_transform=None,
    scaler=None,
    scaler_params=None,
    selector=None,
    selector_params=None,
    classifier=None,
    classifier_params=None):
    from sklearn.pipeline import Pipeline
    import numpy as np
    from estimators import get_splitter, get_classifier, get_selector, get_scaler

    steps = []
    if zero_fraction_filter is not None:
        steps.append(('zero_fraction_filter', get_selector('zero_fraction_filter', 
            **parse_params(zero_fraction_filter))))
    if rpkm_filter != 'none':
        if feature_names is None:
            raise ValueError('feature_names is required for rpkm_filter')
        gene_lengths = get_gene_lengths_from_feature_names(feature_names)
        step = get_selector('rpkm_filter', **parse_params(rpkm_filter))
        step.set_gene_lengths(gene_lengths)
        steps.append(('rpkm_filter', step))
    if fold_change_filter != 'none':
        steps.append(('fold_change_filter', get_selector(
            'fold_change_filter', **parse_params(fold_change_filter))))
    if log_transform != 'none':
        steps.append(('log_transform', get_scaler(
            'log_transform', **parse_params(log_transform))))
    if scaler:
        steps.append(('scaler', get_scaler(scaler, **parse_params(scaler))))
    if selector:
        if selector in ('robust', 'rfe', 'rfe_cv'):
            estimator = get_classifier(classifier, **parse_params(classifier_params))
        else:
            estimator = None
        steps.append(('selector', get_selector(
            selector, estimator=estimator, **parse_params(selector_params))))
    steps.append(('classifier', get_classifier(
        classifier, **parse_params(classifier_params))))
    
    pipeline = Pipeline(steps)
    return pipeline

def read_data_matrix(matrix, sample_classes, transpose=False, positive_class=None, negative_class=None):
    # read data matrix
    logger.info('read data matrix: ' + matrix)
    X = pd.read_table(matrix, index_col=0, sep='\t')
    # transpose
    if transpose:
        logger.info('transpose feature matrix')
        X = X.T
    logger.info('number of features: {}'.format(X.shape[1]))
    # read sample classes
    logger.info('read sample classes: ' + sample_classes)
    sample_classes = pd.read_table(sample_classes, index_col=0, sep='\t')
    sample_classes = sample_classes.iloc[:, 0]
    sample_classes = sample_classes.loc[X.index.values]
    logger.info('sample_classes: {}'.format(sample_classes.shape[0]))
    
    # get positive and negative classes
    if (positive_class is not None) and (negative_class is not None):
        positive_class = positive_class.split(',')
        negative_class = negative_class.split(',')
    else:
        unique_classes = np.unique(sample_classes.values)
        if len(unique_classes) != 2:
            raise ValueError('expect 2 classes but {} classes found'.format(len(unique_classes)))
        positive_class, negative_class = unique_classes
    positive_class = np.atleast_1d(positive_class)
    negative_class = np.atleast_1d(negative_class)
    # select positive samples and negative samples
    logger.info('positive class: {}'.format(positive_class))
    logger.info('negative class: {}'.format(negative_class))
    X_pos = X.loc[sample_classes[sample_classes.isin(positive_class)].index.values]
    X_neg = X.loc[sample_classes[sample_classes.isin(negative_class)].index.values]
    logger.info('number of positive samples: {}, negative samples: {}, class ratio: {}'.format(
        X_pos.shape[0], X_neg.shape[0], float(X_pos.shape[0])/X_neg.shape[0]))
    X = pd.concat([X_pos, X_neg], axis=0)
    y = np.zeros(X.shape[0], dtype=np.int32)
    y[X_pos.shape[0]:] = 1
    del X_pos
    del X_neg
    n_samples, n_features = X.shape
    sample_ids = X.index.values
    feature_names = X.columns.values
    X = X.values

    return X, y, sample_ids, feature_names

@command_handler
def cross_validation(args):
    from estimators2 import search_dict, CollectMetrics, CollectPredictions, FeatureSelectionMatrix,\
        CombinedEstimator, get_features_from_pipeline
    from estimators2 import cross_validation as _cross_validation
    import pandas as pd
    import numpy as np

    logger.info('parameters: {}'.format(vars(args)))
    # read data matrix
    X, y, sample_ids, feature_names = read_data_matrix(args.matrix, args.sample_classes,
        **search_dict(vars(args), ('transpose', 'positive_class', 'negative_class')))
    
    #estimator = make_estimator(feature_names,
    #    **search_dict(vars(args), ('zero_fraction_filter', 'rpkm_filter', 'fold_change_filter',
    #    'log_transform', 'scaler', 'scaler_params',
    #    'selector', 'selector_params', 'classifier', 'classifier_params')))
    estimator = CombinedEstimator(feature_names,
        **search_dict(vars(args), ('zero_fraction_filter', 'rpkm_filter', 'fold_change_filter',
        'log_transform', 'scaler', 'scaler_params',
        'selector', 'selector_params', 'classifier', 'classifier_params',
        'grid_search', 'grid_search_param_grid', 'grid_search_cv_params', 'grid_search_scoring')))
    
    #for name in estimator.named_steps:
    #    logger.info('step: {}'.format(name))

    logger.info('start cross-validation')
    cv_callbacks = [CollectMetrics(), CollectPredictions(), FeatureSelectionMatrix()]
    cv_params = parse_params(args.cv_params)
    _cross_validation(estimator, X, y, params=cv_params, callbacks=cv_callbacks)
    logger.info('collect_metrics:')
    print(cv_callbacks[0].get_metrics())
    logger.info('fit estimator on full dataset')
    estimator.fit(X, y)
    logger.info('classifier params: {}'.format(estimator.classifier_.get_params()))
    if args.selector is not None:
        feature_index = estimator.features_
        logger.info('number of selected features: {}'.format(feature_index.shape[0]))
        logger.info('selected features')
        print(feature_names[feature_index])

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Machine learning module')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('cross_validation')
    g_input = parser.add_argument_group('input')
    g_input.add_argument('--matrix', '-i', type=str, metavar='FILE', required=True,
        help='input feature matrix (rows are samples and columns are features')
    g_input.add_argument('--sample-classes', type=str, metavar='FILE', required=True,
        help='input file containing sample classes with 2 columns: sample_id, sample_class')
    g_input.add_argument('--positive-class', type=str, metavar='STRING',
        help='comma-separated list of sample classes to use as positive class')
    g_input.add_argument('--negative-class', type=str,metavar='STRING',
        help='comma-separates list of sample classes to use as negative class')
    g_input.add_argument('--transpose', action='store_true', default=False,
        help='transpose the feature matrix')

    g_filter = parser.add_argument_group('filter')
    g_filter.add_argument('--zero-fraction-filter', type=str, nargs='?', default='none')
    g_filter.add_argument('--rpkm-filter', type=str, nargs='?', default='none')
    g_filter.add_argument('--fold-change-filter', type=str, nargs='?', default='none')

    g_scaler = parser.add_argument_group('scaler')
    g_scaler.add_argument('--log-transform', type=str, nargs='?', default='none')
    g_scaler.add_argument('--scaler', type=str, metavar='NAME', default='robust')
    g_scaler.add_argument('--scaler-params', type=str, metavar='STRING')

    g_select = parser.add_argument_group('feature_selection')
    g_select.add_argument('--selector', type=str, metavar='NAME', default='robust')
    g_select.add_argument('--selector-params', type=str, metavar='STRING')

    g_classifier = parser.add_argument_group('classifier')
    g_classifier.add_argument('--classifier', type=str, metavar='NAME', default='random_forest')
    g_classifier.add_argument('--classifier-params', type=str, metavar='STRING')

    g_cv = parser.add_argument_group('cross_validation')
    g_cv.add_argument('--cv-params', type=str, metavar='STRING', nargs='?')
    g_cv.add_argument('--grid-search', action='store_true')
    g_cv.add_argument('--grid-search-param-grid', type=str, metavar='STRING')
    g_cv.add_argument('--grid-search-cv-params', type=str, metavar='STRING')
    g_cv.add_argument('--grid-search-scoring', type=str, metavar='STRING')

    g_misc= parser.add_argument_group('misc')
    g_misc.add_argument('--compute-sample-weight', action='store_true',
        help='compute sample weight to balance classes')
    
    g_output= parser.add_argument_group('output')
    g_output.add_argument('--output-dir', '-o', type=str, metavar='DIR', 
        required=True, help='output directory')
    
    args = main_parser.parse_args()
    if args.command is None:
        print('Errror: missing command', file=sys.stdout)
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('machine_learning.' + args.command)

    import pandas as pd
    import numpy as np

    command_handlers.get(args.command)(args)