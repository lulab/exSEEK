#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

@command_handler
def evaluate(args):
    import numpy as np
    import pandas as pd
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import LinearSVC
    from sklearn.metrics import roc_auc_score, accuracy_score, get_scorer
    from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler, MaxAbsScaler
    from sklearn.model_selection import GridSearchCV
    from sklearn.feature_selection import RFE, RFECV
    from sklearn.utils.class_weight import compute_sample_weight
    from sklearn.model_selection import KFold, StratifiedKFold, ShuffleSplit, LeaveOneOut, \
        RepeatedKFold, RepeatedStratifiedKFold, LeaveOneOut, StratifiedShuffleSplit
    import pickle
    from estimators import RobustEstimator
    from tqdm import tqdm
    import h5py

    logger.info('read feature matrix: ' + args.matrix)
    m = pd.read_table(args.matrix, index_col=0, sep='\t')
    if args.transpose:
        logger.info('transpose feature matrix')
        m = m.T
    logger.info('{} samples, {} features'.format(m.shape[0], m.shape[1]))
    if args.remove_zero_features is not None:
        logger.info('remove features with zero fraction below {}'.format(args.remove_zero_features))
        m = m.loc[:, np.isclose(m, 0).sum(axis=0) < (m.shape[0]*args.remove_zero_features)]
    if args.top_features_by_median is not None:
        nonzero_samples = (m > 0).sum(axis=0)
        counts_geomean = np.exp(np.sum(np.log(np.maximum(m, 1)), axis=0)/nonzero_samples)
        m = m.loc[:, counts_geomean.sort_values(ascending=False)[:args.top_features_by_median].index.values]
    feature_names = m.columns.values
    logger.info('{} samples, {} features'.format(m.shape[0], m.shape[1]))
    logger.info('sample: {} ...'.format(str(m.index.values[:3])))
    logger.info('features: {} ...'.format(str(m.columns.values[:3])))

    logger.info('read sample classes: ' + args.sample_classes)
    sample_classes = pd.read_table(args.sample_classes, header=None, 
        names=['sample_id', 'sample_class'], index_col=0, sep='\t')
    sample_classes = sample_classes.iloc[:, 0]
    sample_classes = sample_classes.loc[m.index.values]
    logger.info('sample_classes: {}'.format(sample_classes.shape[0]))

    # select samples
    if (args.positive_class is not None) and (args.negative_class is not None):
        positive_class = args.positive_class.split(',')
        negative_class = args.negative_class.split(',')
    else:
        unique_classes = np.unique(sample_classes.values)
        if len(unique_classes) != 2:
            raise ValueError('expect 2 classes but {} classes found'.format(len(unique_classes)))
        positive_class, negative_class = unique_classes
    positive_class = np.atleast_1d(positive_class)
    negative_class = np.atleast_1d(negative_class)

    logger.info('positive class: {}, negative class: {}'.format(positive_class, negative_class))
    X_pos = m.loc[sample_classes[sample_classes.isin(positive_class)].index.values]
    X_neg = m.loc[sample_classes[sample_classes.isin(negative_class)].index.values]
    logger.info('number of positive samples: {}, negative samples: {}, class ratio: {}'.format(
        X_pos.shape[0], X_neg.shape[0], float(X_pos.shape[0])/X_neg.shape[0]))
    X = pd.concat([X_pos, X_neg], axis=0)
    y = np.zeros(X.shape[0], dtype=np.int32)
    y[X_pos.shape[0]:] = 1
    del X_pos
    del X_neg
    X_raw = X
    n_samples, n_features = X.shape
    sample_ids = X.index.values

    if not os.path.isdir(args.output_dir):
        logger.info('create outout directory: ' + args.output_dir)
        os.makedirs(args.output_dir)

    logger.info('save sample ids')
    X.index.to_series().to_csv(os.path.join(args.output_dir, 'samples.txt'),
        sep='\t', header=False, index=False)
    logger.info('save sample classes')
    np.savetxt(os.path.join(args.output_dir, 'classes.txt'), y, fmt='%d')

    if args.use_log:
        logger.info('apply log2 to feature matrix')
        X = np.log2(X + 0.001)

    if args.scaler == 'zscore':
        logger.info('scale features using z-score normalization')
        X = StandardScaler().fit_transform(X)
    elif args.scaler == 'robust':
        logger.info('scale features using robust normalization')
        X = RobustScaler().fit_transform(X)
    elif args.scaler == 'min_max':
        logger.info('scale features using min-max normalization')
        X = MinMaxScaler().fit_transform(X)
    elif args.scaler == 'max_abs':
        logger.info('scale features using max-abs normalization')
        X = MaxAbsScaler().fit_transform(X)

    # check NaN values
    if np.any(np.isnan(X)):
        logger.info('nan values found in features')
    estimator = None
    grid_search = None
    logger.info('use {} to select features'.format(args.method))
    if args.method == 'logistic_regression':
        estimator = LogisticRegression()
        grid_search = {'C': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2, 1e3, 1e4, 1e5]}
    elif args.method == 'random_forest':
        estimator = RandomForestClassifier()
        grid_search = {'n_estimators': [10, 20, 30], 'max_depth': list(range(2, 10)) }
    elif args.method == 'linear_svm':
        estimator = LinearSVC()
        grid_search = {'C': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2, 1e3, 1e4, 1e5]}
    else:
        raise ValueError('unknown feature selection method: {}'.format(args.method))
    
    #logger.info('load model: ' + args.model_file)
    #with open(args.model_file, 'rb') as f:
    #    model = pickle.load(f)
    
    def get_splitter(splitter, n_splits=5, n_repeats=5, test_size=0.2):
        if splitter == 'kfold':
            return KFold(n_splits=n_splits)
        elif splitter == 'stratified_kfold':
            return StratifiedKFold(n_splits=n_splits)
        elif splitter == 'repeated_stratified_kfold':
            return RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats)
        elif splitter == 'shuffle_split':
            return ShuffleSplit(n_splits=n_splits, test_size=test_size)
        elif splitter == 'stratified_shuffle_split':
            return StratifiedShuffleSplit(n_splits=n_splits, test_size=test_size)
        elif splitter == 'leave_one_out':
            return LeaveOneOut()
        else:
            raise ValueError('unknown splitter: {}'.format(splitter))

    def score_function(estimator):
        '''Get method of an estimator that predict a continous score for each sample
        '''
        if hasattr(estimator, 'predict_proba'):
            return estimator.predict_proba
        elif hasattr(estimator, 'decision_function'):
            return estimator.decision_function
        else:
            raise ValueError('the estimator should either have decision_function() method or predict_proba() method')

    def feature_importances(estimator):
        '''Get feature importance attribute of an estimator
        '''
        if hasattr(estimator, 'coef_'):
            return np.ravel(estimator.coef_)
        elif hasattr(estimator, 'feature_importances_'):
            return np.ravel(estimator.feature_importances_)
        else:
            raise ValueError('the estimator should have either coef_ or feature_importances_ attribute')
    
    def get_scorer(scoring):
        if scoring == 'roc_auc':
            return roc_auc_score
        else:
            raise ValueError('unknonwn scoring: {}'.format(scoring))

    splitter = get_splitter(args.splitter, n_splits=args.n_splits, n_repeats=args.n_repeats)
    metrics = []
    predictions = np.full((splitter.get_n_splits(X), X.shape[0]), np.nan)
    predicted_labels = np.full((splitter.get_n_splits(X), X.shape[0]), np.nan)
    train_index_matrix = np.zeros((splitter.get_n_splits(X), X.shape[0]),dtype=np.bool)
    feature_selection_matrix = None
    if args.n_select is not None:
        feature_selection_matrix = np.zeros((splitter.get_n_splits(X), X.shape[1]), dtype=bool)
    if args.rfe:
        if 0.0 < args.rfe_step < 1.0:
            rfe_step = int(max(1, args.rfe_step*n_features))
        else:
            rfe_step = int(args.rfe_step)
        rfe_scores = None
    i_split = 0
    scorer = get_scorer(args.scorer)
    data_splits = list(splitter.split(X, y))
    data_splits.append((np.arange(n_samples), None))
    for train_index, test_index in tqdm(data_splits, total=splitter.get_n_splits(X) + 1, unit='fold'):
        X_train, y_train = X[train_index], y[train_index]
        X_test, y_test = X[test_index], y[test_index]
        # optimize hyper-parameters
        if grid_search is not None:
            cv = GridSearchCV(estimator, grid_search, cv=3)
            cv.fit(X[train_index], y[train_index])
            estimator = cv.best_estimator_
        
        sample_weight = np.ones(X_train.shape[0])
        if args.compute_sample_weight:
            sample_weight = compute_sample_weight('balanced', y_train)
        # feature selection
        if args.n_select is not None:
            if args.robust_select:
                resampler_args = {}
                if args.robust_resample_method == 'jackknife':
                    resampler_args = {'max_runs': args.robust_max_runs,
                        'remove': args.robust_jackknife_remove
                    }
                elif args.robust_resample_method == 'bootstrap':
                    resampler_args = {'max_runs': args.robust_max_runs}
                robust_estimator = RobustEstimator(estimator, n_select=args.n_select, 
                    grid_search=grid_search, resample_method=args.robust_resample_method,
                    rfe=args.rfe, **resampler_args)
                robust_estimator.fit(X_train, y_train, sample_weight=sample_weight)
                estimator = robust_estimator.estimator_
                features = robust_estimator.features_
            # RFE feature selection
            elif args.rfe:
                rfe = RFE(estimator, n_features_to_select=args.n_select, step=rfe_step)
                if i_split < splitter.get_n_splits(X):
                    if args.splitter == 'leave_one_out':
                        # AUC is undefined for only one test sample
                        step_score = lambda estimator, features: np.nan
                    else:
                        step_score = lambda estimator, features: scorer(y_test, 
                            score_function(estimator)(X[test_index][:, features])[:, 1])
                else:
                    step_score = None
                rfe._fit(X_train, y_train, step_score=step_score)
                features = np.nonzero(rfe.ranking_ == 1)[0]
                if i_split < splitter.get_n_splits(X):
                    if rfe_scores is None:
                        rfe_n_steps = len(rfe.scores_)
                        rfe_n_features_step = np.maximum(n_features - rfe_step*np.arange(rfe_n_steps), 1)
                        rfe_scores = np.zeros((splitter.get_n_splits(X), rfe_n_steps))
                    rfe_scores[i_split] = rfe.scores_
                estimator = rfe.estimator_
            else:
                # train the model
                estimator.fit(X[train_index], y[train_index], sample_weight=sample_weight)
                features = np.argsort(-feature_importances(estimator))[:args.n_select]
            if i_split < splitter.get_n_splits(X):
                feature_selection_matrix[i_split, features] = True
        else:
            # no feature selection
            features = np.arange(n_features, dtype=np.int64)
        
        estimator.fit(X[train_index][:, features], y[train_index], sample_weight=sample_weight)
        if i_split != splitter.get_n_splits(X):
            predictions[i_split] = score_function(estimator)(X[:, features])[:, 1]
            predicted_labels[i_split] = estimator.predict(X[:, features])
            metric = {}
            metric['train_{}'.format(args.scorer)] = scorer(y_train, predictions[i_split, train_index])
            # AUC is undefined for only one test sample
            if args.splitter != 'leave_one_out':
                metric['test_{}'.format(args.scorer)] = scorer(y_test, predictions[i_split, test_index])
            if args.splitter in ('repeated_kfold', 'repeated_stratified_kfold'):
                metric['repeat'] = i_split//args.n_repeats
                metric['split'] = i_split%args.n_repeats
            else:
                metric['split'] = i_split
            metrics.append(metric)
            train_index_matrix[i_split, train_index] = True
        i_split += 1
    metrics = pd.DataFrame.from_records(metrics)
    if args.splitter == 'leave_one_out':
        metrics['test_{}'.format(args.scorer)] = scorer(y, predictions[np.r_[:n_samples], np.r_[:n_samples]])

    logger.info('save best model')
    with open(os.path.join(args.output_dir, 'best_model.pkl'), 'wb') as f:
        pickle.dump(estimator, f)
    
    logger.info('save features')
    data = pd.Series(features, index=feature_names[features])
    data.to_csv(os.path.join(args.output_dir, 'features.txt'), sep='\t', header=False)

    logger.info('save feature importances')
    data = pd.Series(feature_importances(estimator), index=feature_names[features])
    data.to_csv(os.path.join(args.output_dir, 'feature_importances.txt'), sep='\t', header=False)

    logger.info('save evaluations')
    with h5py.File(os.path.join(args.output_dir, 'evaluation.{}.h5'.format(args.splitter)), 'w') as f:
        f.create_dataset('train_index', data=train_index_matrix)
        f.create_dataset('predictions', data=predictions)
        if feature_selection_matrix is not None:
            f.create_dataset('feature_selection', data=feature_selection_matrix)
        if args.rfe:
            f.create_dataset('rfe_n_features_step', data=rfe_n_features_step)
            f.create_dataset('rfe_scores', data=rfe_scores)
        f.create_dataset('labels', data=y)
        f.create_dataset('predicted_labels', data=predicted_labels)
        
    logger.info('save metrics')
    metrics.to_csv(os.path.join(args.output_dir, 'metrics.{}.txt'.format(args.splitter)), sep='\t', header=True, index=False)

@command_handler
def robust_select(args):
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import seaborn as sns
    sns.set()
    from sklearn.linear_model import LogisticRegression
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.svm import LinearSVC
    from sklearn.metrics import roc_auc_score, accuracy_score
    from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler, MaxAbsScaler
    from sklearn.utils.class_weight import compute_sample_weight
    from sklearn.model_selection import GridSearchCV
    import os
    from tqdm import tqdm
    import pickle
    import json
    from scipy.cluster import hierarchy
    from scipy.spatial.distance import pdist, squareform
    from estimators import RobustEstimator, TTestEstimator

    logger.info('read feature matrix: ' + args.matrix)
    m = pd.read_table(args.matrix, index_col=0)
    if args.transpose:
        logger.info('transpose feature matrix')
        m = m.T
    if args.remove_zero_features is not None:
        m = m.loc[:, np.isclose(m, 0).sum(axis=0) >= (m.shape[0]*args.remove_zero_features)]
    if args.top_features_by_median is not None:
        nonzero_samples = (m > 0).sum(axis=0)
        counts_geomean = np.exp(np.sum(np.log(np.maximum(m, 1)), axis=0)/nonzero_samples)
        m = m.loc[:, counts_geomean.sort_values(ascending=False)[:args.top_features_by_median].index.values]
    feature_names = m.columns.values
    logger.info('{} samples, {} features'.format(m.shape[0], m.shape[1]))
    logger.info('sample: {} ...'.format(str(m.index.values[:3])))
    logger.info('features: {} ...'.format(str(m.columns.values[:3])))

    logger.info('read sample classes: ' + args.sample_classes)
    sample_classes = pd.read_table(args.sample_classes, header=None, 
        names=['sample_id', 'sample_class'], index_col=0)
    sample_classes = sample_classes.iloc[:, 0]
    sample_classes = sample_classes.loc[m.index.values]
    print('sample_classes: {}'.format(sample_classes.shape[0]))

    # select samples
    if (args.positive_class is not None) and (args.negative_class is not None):
        positive_class = args.positive_class
        negative_class = args.negative_class
    else:
        unique_classes = np.unique(sample_classes.values)
        if len(unique_classes) != 2:
            raise ValueError('expect 2 classes but {} classes found'.format(len(unique_classes)))
        positive_class, negative_class = unique_classes
    logger.info('positive class: {}, negative class: {}'.format(positive_class, negative_class))
    X_pos = m.loc[sample_classes[sample_classes == positive_class].index.values]
    X_neg = m.loc[sample_classes[sample_classes == negative_class].index.values]
    logger.info('number of positive samples: {}, negative samples: {}, class ratio: {}'.format(
        X_pos.shape[0], X_neg.shape[0], float(X_pos.shape[0])/X_neg.shape[0]))
    X = pd.concat([X_pos, X_neg], axis=0)
    y = np.zeros(X.shape[0], dtype=np.int32)
    y[X_pos.shape[0]:] = 1
    del X_pos
    del X_neg
    X_raw = X

    if not os.path.isdir(args.output_dir):
        logger.info('create outout directory: ' + args.output_dir)
        os.makedirs(args.output_dir)

    if args.use_log:
        logger.info('apply log2 to feature matrix')
        X = np.log2(X + 1)

    if args.scaler == 'zscore':
        logger.info('scale features using z-score normalization')
        X = StandardScaler().fit_transform(X)
    elif args.scaler == 'robust':
        logger.info('scale features using robust normalization')
        X = RobustScaler().fit_transform(X)
    elif args.scaler == 'min_max':
        logger.info('scale features using min-max normalization')
        X = MinMaxScaler().fit_transform(X)
    elif args.scaler == 'max_abs':
        logger.info('scale features using max-abs normalization')
        X = MaxAbsScaler().fit_transform(X)

    if np.any(np.isnan(X)):
        logger.info('nan values found in features')
    estimator = None
    grid_search = None
    logger.info('use {} to select features'.format(args.method))
    if args.method == 'r_test':
        estimator = TTestEstimator()
    elif args.method == 'logistic_regression':
        estimator = LogisticRegression()
        grid_search = {'C': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2, 1e3, 1e4, 1e5]}
    elif args.method == 'random_forest':
        estimator = RandomForestClassifier()
        grid_search = {'max_depth': list(range(2, 10))}
    elif args.method == 'linear_svm':
        estimator = LinearSVC()
        grid_search = {'C': [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2, 1e3, 1e4, 1e5]}
    else:
        raise ValueError('unknown feature selection method: {}'.format(args.method))

    resampler_args = {}
    if args.resample_method == 'jackknife':
        resampler_args = {'max_runs': args.jackknife_max_runs,
            'remove': args.jackknife_remove
        }
    elif args.resample_method == 'bootstrap':
        resampler_args = {'max_runs': args.bootstrap_max_runs}

    robust_estimator = RobustEstimator(estimator, n_select=args.n_select, 
        grid_search=grid_search, resample_method='bootstrap',
        rfe=args.rfe, **resampler_args)

    sample_weight = None
    if args.compute_sample_weight:
        sample_weight = compute_sample_weight('balanced', y)
    robust_estimator.fit(X, y, sample_weight=sample_weight)

    logger.info('save model')
    with open(os.path.join(args.output_dir, 'model.pkl'), 'wb') as f:
        pickle.dump(robust_estimator.estimator_, f)
    
    logger.info('save model parameters')
    with open(os.path.join(args.output_dir, 'params.json'), 'w') as f:
        json.dump(robust_estimator.estimator_.get_params(), f, indent=2)

    logger.info('plot histogram for feature recurrence')
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.hist(robust_estimator.feature_recurrence_, bins=30)
    ax.set_xlabel('Occurence')
    ax.set_ylabel('Counts')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'feature_recurrence_hist.pdf'))
    plt.close()

    """
    logger.info('plot heatmap of raw features')
    sns.clustermap(np.log2(X_raw + 1), cmap='vlag', figsize=(16, 20))
    plt.savefig(os.path.join(args.output_dir, 'heatmap_raw.pdf'))
    plt.close()

    logger.info('plot heatmap scaled features')
    sns.clustermap(pd.DataFrame(X, index=X_raw.index, columns=X_raw.columns), cmap='vlag', figsize=(16, 20))
    plt.savefig(os.path.join(args.output_dir, 'heatmap_scaled.pdf'))
    plt.close()
    """

    logger.info('plot recurrence of robust features')
    df = pd.DataFrame({'Feature': feature_names[robust_estimator.features_],
                   'Feature recurrence': robust_estimator.feature_recurrence_[robust_estimator.features_]})
    fig, ax = plt.subplots(figsize=(14, 10))
    sns.barplot('Feature recurrence', 'Feature', color='gray',
            data=df, ax=ax)
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'feature_recurrence_barplot.pdf'))
    plt.close()

    logger.info('plot heatmap of feature recurrence')
    fig, ax = plt.subplots(figsize=(14, 14))
    cmap = sns.light_palette((127/360.0, 40/100.0, 49/100.0), n_colors=2, input='hls', as_cmap=True)
    logger.info('number of features: {}'.format(len(robust_estimator.features_)))
    df = pd.DataFrame(robust_estimator.feature_selection_matrix_,
        index=np.arange(1, robust_estimator.feature_selection_matrix_.shape[0] + 1), 
        columns=feature_names[robust_estimator.features_]).T
    sns.heatmap(df, linewidth=0.1, cmap=cmap, ax=ax, cbar=False)
    ax.set_xlabel('Run')
    ax.set_ylabel('Feature')
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, 'feature_recurrence_heatmap.pdf'))
    plt.close()

    logger.info('save feature selection matrix')
    data = pd.DataFrame(robust_estimator.feature_selection_matrix_,
        index=np.arange(1, robust_estimator.feature_selection_matrix_.shape[0] + 1),
        columns=feature_names[robust_estimator.features_]).T
    data.index.name = 'run'
    data.columns.name = 'peak'
    data.to_csv(os.path.join(args.output_dir, 'feature_selection_matrix.txt'), sep='\t', header=True, index=True)

    logger.info('plot feature importances of refitted model')
    fig, ax = plt.subplots(figsize=(12, 10))
    df = pd.DataFrame({'feature': feature_names[robust_estimator.features_],
                    'feature_importance': robust_estimator.feature_importances_})
    sns.barplot('feature_importance', 'feature', color='gray',
            data=df, ax=ax)
    plt.savefig(os.path.join(args.output_dir, 'feature_importances_refitted.pdf'))
    plt.close()

    logger.info('plot heatmap of selected features (scaled)')
    df = pd.DataFrame(X, index=X_raw.index, columns=feature_names)
    df = df.iloc[:, robust_estimator.features_].T
    sns.clustermap(df, cmap='vlag', figsize=(15, 12), row_cluster=True)
    plt.subplots_adjust(right=0.8, bottom=0.2)
    plt.savefig(os.path.join(args.output_dir, 'heatmap_selected_features_scaled.pdf'))
    plt.close()

    logger.info('plot correlation between selected features (scaled)')
    df = pd.DataFrame(X, index=X_raw.index, columns=feature_names)
    df = df.iloc[:, robust_estimator.features_].T
    distmat = np.nan_to_num(1 - squareform(pdist(df.values, metric='correlation')))
    distmat = pd.DataFrame(distmat, index=df.index, columns=df.index)
    sns.clustermap(distmat, cmap='vlag', figsize=(15, 15), row_cluster=True)
    plt.subplots_adjust(right=0.8, bottom=0.2)
    plt.savefig(os.path.join(args.output_dir, 'heatmap_selected_feature_correlation.pdf'))
    plt.close()

    logger.info('plot correlation between samples (scaled)')
    df = pd.DataFrame(X, index=X_raw.index, columns=feature_names)
    df = df.iloc[:, robust_estimator.features_]
    distmat = np.nan_to_num(1 - squareform(pdist(df.values, metric='correlation')))
    distmat = pd.DataFrame(distmat, index=df.index, columns=df.index)
    sns.clustermap(distmat, cmap='vlag', figsize=(15, 15), row_cluster=True)
    plt.subplots_adjust(right=0.8, bottom=0.2)
    plt.savefig(os.path.join(args.output_dir, 'heatmap_sample_correlation.pdf'))
    plt.close()

    logger.info('plot heatmap of selected features (raw)')
    df = X_raw.iloc[:, robust_estimator.features_].T
    g = sns.clustermap(np.log2(df + 1), cmap='vlag', figsize=(15, 12), row_cluster=True)
    plt.subplots_adjust(right=0.8, bottom=0.2)
    plt.savefig(os.path.join(args.output_dir, 'heatmap_selected_features_raw.pdf'))
    plt.close()

    logger.info('save feature importance matrix')
    data = pd.DataFrame(robust_estimator.feature_importances_matrix_, columns=feature_names)
    data.to_csv(os.path.join(args.output_dir, 'feature_importance_matrix.txt'), sep='\t', index=False, header=True, na_rep='nan')

    logger.info('save feature recurrence matrix')
    data = pd.Series(robust_estimator.feature_recurrence_, index=feature_names)
    data.index.name = 'peak'
    data.to_csv(os.path.join(args.output_dir, 'feature_recurrence.txt'), sep='\t', index=True, header=False, na_rep='nan')

    logger.info('save feature importances estimates (mean and std)')
    data = pd.DataFrame({'feature_importances_mean': robust_estimator.feature_importances_mean_,
         'feature_importances_std': robust_estimator.feature_importances_std_},
         index=feature_names)
    data.index.name = 'peak'
    data.to_csv(os.path.join(args.output_dir, 'feature_importances_estimate.txt'), sep='\t', index=True, header=True, na_rep='nan')

    logger.info('save selected features')
    data = pd.Series(robust_estimator.features_, index=feature_names[robust_estimator.features_])
    data.to_csv(os.path.join(args.output_dir, 'features.txt'), header=False, index=True, sep='\t')

    logger.info('save selected feature importances')
    data = pd.Series(robust_estimator.feature_importances_, index=feature_names[robust_estimator.features_])
    data.to_csv(os.path.join(args.output_dir, 'feature_importances.txt'), header=False, index=True, sep='\t')

@command_handler
def preprocess(args):
    from estimators import FeatureSelection

    fs = FeatureSelection(**vars(args))
    fs.read_data()
    fs.filter_features()
    fs.select_samples()
    fs.scale_features()
    fs.save_matrix()
    fs.single_feature_metrics()
    fs.plot_feature_importance()
    fs.load_model()
    fs.evaluate_model()

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Feature selection module')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('robust_select')
    parser.add_argument('--matrix', '-i', type=str, required=True,
        help='input feature matrix (rows are samples and columns are features')
    parser.add_argument('--sample-classes', type=str, required=True,
        help='input file containing sample classes with 2 columns: sample_id, sample_class')
    parser.add_argument('--remove-zero-features', type=float, 
        help='remove features that have fraction of zero values above this value')
    parser.add_argument('--top-features-by-median', type=int,
        help='select this number of features with highest median values across samples')
    parser.add_argument('--transpose', action='store_true', default=False,
        help='transpose the feature matrix')
    parser.add_argument('--positive-class', type=str,
        help='sample class to use as positive class')
    parser.add_argument('--negative-class', type=str,
        help='sample class to use as negative class')
    parser.add_argument('--use-log', action='store_true',
        help='apply log2 to feature matrix')
    parser.add_argument('--scaler', type=str, default='zscore',
        choices=('zscore', 'robust', 'max_abs', 'min_max', 'none'),
        help='method for scaling features')
    parser.add_argument('--method', type=str, default='logistic_regression',
        choices=('t_test', 'logistic_regression', 'random_forest', 'linear_svm'),
        help='feature selection method')
    parser.add_argument('--n-select', type=int,
        help='number of features to select')
    parser.add_argument('--output-dir', '-o', type=str, required=True,
        help='output directory')
    parser.add_argument('--compute-sample-weight', action='store_true',
        help='compute sample weight to balance classes')
    parser.add_argument('--rfe', action='store_true', default=False,
        help='use RFE to select features')
    parser.add_argument('--resample-method', type=str, default='jackknife',
        choices=('jackknife', 'bootstrap'), help='resampling method')
    parser.add_argument('--jackknife-max-runs', type=int, default=100)
    parser.add_argument('--jackknife-remove', type=int, default=0.2)
    parser.add_argument('--bootstrap-max-runs', type=int, default=100)

    parser = subparsers.add_parser('preprocess')
    parser.add_argument('--matrix', '-i', type=str, required=True,
        help='input feature matrix (rows are samples and columns are features')
    parser.add_argument('--sample-classes', type=str, required=True,
        help='input file containing sample classes with 2 columns: sample_id, sample_class')
    parser.add_argument('--remove-zero-features', type=float, 
        help='remove features that have fraction of zero values above this value')
    parser.add_argument('--top-features-by-median', type=int,
        help='select this number of features with highest median values across samples')
    parser.add_argument('--transpose', action='store_true', default=False,
        help='transpose the feature matrix')
    parser.add_argument('--positive-class', type=str,
        help='comma-separated list of sample classes to use as positive class')
    parser.add_argument('--negative-class', type=str,
        help='comma-separates list of sample classes to use as negative class')
    parser.add_argument('--use-log', action='store_true',
        help='apply log2 to feature matrix')
    parser.add_argument('--scaler', type=str, default='zscore',
        choices=('zscore', 'robust', 'max_abs', 'min_max', 'none'),
        help='method for scaling features')
    parser.add_argument('--output-dir', '-o', type=str, required=True,
        help='output directory')

    parser = subparsers.add_parser('evaluate')
    parser.add_argument('--matrix', '-i', type=str, required=True,
        help='input feature matrix (rows are samples and columns are features')
    parser.add_argument('--sample-classes', type=str, required=True,
        help='input file containing sample classes with 2 columns: sample_id, sample_class')
    parser.add_argument('--transpose', action='store_true', default=False,
        help='transpose the feature matrix')
    parser.add_argument('--positive-class', type=str,
        help='comma-separated list of sample classes to use as positive class')
    parser.add_argument('--negative-class', type=str,
        help='comma-separates list of sample classes to use as negative class')
    parser.add_argument('--use-log', action='store_true',
        help='apply log2 to feature matrix')
    parser.add_argument('--scaler', type=str, default='zscore',
        choices=('zscore', 'robust', 'max_abs', 'min_max', 'none'),
        help='method for scaling features')
    parser.add_argument('--method', type=str, default='logistic_regression',
        choices=('logistic_regression', 'random_forest', 'linear_svm'),
        help='feature selection method')
    parser.add_argument('--output-dir', '-o', type=str, required=True,
        help='output directory')
    parser.add_argument('--compute-sample-weight', action='store_true',
        help='compute sample weight to balance classes')
    parser.add_argument('--splitter', type=str, default='stratified_kfold',
        choices=('kfold', 'stratified_kfold', 'shuffle_split', 'repeated_stratified_kfold', 'leave_one_out', 'stratified_shuffle_split'))
    parser.add_argument('--rfe', action='store_true', default=False,
        help='use RFE to select features')
    parser.add_argument('--rfe-step', type=float,
        help='number/fraction of features to eliminate in each step')
    parser.add_argument('--n-select', type=int,
        help='number of features to select')
    parser.add_argument('--n-splits', type=int, default=5,
        help='number of splits for kfold, stratified_kfold and shuffle_splits')
    parser.add_argument('--n-repeats', type=int, default=10,
        help='number of repeats for repeated_stratified_kfold and repeated_kfold')
    parser.add_argument('--test-size', type=float, default=0.2,
        help='fraction/number of samples for testing')
    parser.add_argument('--scorer', type=str, default='roc_auc',
        help='metric to use')
    parser.add_argument('--robust-select', action='store_true', default=False,
        help='use robust feature selection by selecting recurring features across resampling runs')
    parser.add_argument('--robust-resample-method', type=str, default='jackknife',
        choices=('bootstrap', 'jackknife'), help='resampling method for robust feature selection')
    parser.add_argument('--robust-max-runs', type=int, default=50,
        help='number of resampling runs for robust feature selections')
    parser.add_argument('--robust-jackknife-remove', type=float, default=1,
        help='number/fraction of samples to remove during each resampling run for robust feature selection')
    parser.add_argument('--top-features-by-median', type=int,
        help='select this number of features with highest median values across samples')
    parser.add_argument('--remove-zero-features', type=float, 
        help='remove features that have fraction of zero values above this value')

    args = main_parser.parse_args()
    if args.command is None:
        print('Errror: missing command', file=sys.stdout)
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('feature_selection.' + args.command)

    command_handlers.get(args.command)(args)