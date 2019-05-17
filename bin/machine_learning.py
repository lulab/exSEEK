#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.DEBUG, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('machine_learning')
import json
from tqdm import tqdm

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

def read_data_matrix(matrix, sample_classes, features=None, transpose=False, positive_class=None, negative_class=None):
    import pandas as pd
    import numpy as np
    # read data matrix
    logger.info('read data matrix: ' + matrix)
    X = pd.read_table(matrix, index_col=0, sep='\t')
    # transpose
    if transpose:
        logger.info('transpose feature matrix')
        X = X.T
    if features is not None:
        logger.info('read subset of feature names from: ' + features)
        features = pd.read_table(features, header=None).iloc[:, 0].values
        logger.info('select {} features'.format(len(features)))
        X = X.reindex(columns=features)
        is_na_features = X.isna().any(axis=0)
        na_features = X.columns.values[is_na_features]
        if na_features.shape[0] > 0:
            logger.warn('missing features found in matrix file: {}'.format(na_features[:10].tolist()))
            #X = X.loc[:, ~is_na_features]
            #raise ValueError('some given features are not found in matrix file: {}'.format(na_features[:10].tolist()))
    logger.info('number of features: {}'.format(X.shape[1]))
    # read sample classes
    logger.info('read sample classes: ' + sample_classes)
    sample_classes = pd.read_table(sample_classes, index_col=0, sep='\t')
    sample_classes = sample_classes.iloc[:, 0]
    sample_classes = sample_classes.loc[X.index.values]
    logger.info('sample_classes: {}'.format(sample_classes.shape[0]))
    
    # get positive and negative classes
    if (positive_class is not None) and (negative_class is not None):
        #positive_class = positive_class.split(',')
        #negative_class = negative_class.split(',')
        pass
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
    # set negative class to 0 and positive class to 1
    y = np.zeros(X.shape[0], dtype=np.int32)
    y[:X_pos.shape[0]] = 1
    del X_pos
    del X_neg
    n_samples, n_features = X.shape
    sample_ids = X.index.values
    feature_names = X.columns.values
    X = X.values

    return X, y, sample_ids, feature_names

def search_params_in_args(args, prefix):
    params = {}
    for key, val in args.items():
        if key.startswith(prefix) and (val is not None):
            params[key[len(prefix):]] = val
    return params

@command_handler
def cross_validation1(args):
    from estimators2 import search_dict, CollectMetrics, CollectPredictions, CollectTrainIndex, FeatureSelectionMatrix,\
        CombinedEstimator, get_features_from_pipeline, parse_params
    from estimators2 import cross_validation as _cross_validation
    import pandas as pd
    import numpy as np
    import h5py
    import pickle
    import yaml

    # read input data matrix
    X, y, sample_ids, feature_names = read_data_matrix(args.matrix, args.sample_classes,
        **search_dict(vars(args), ('features', 'transpose', 'positive_class', 'negative_class')))
    #has_missing_features = np.any(np.isnan(X))
    #if has_missing_features:
        # fill missing features with 0
    #    X[np.isnan(X)] = 0
    # remove missing features
    X = X[:, np.all(~np.isnan(X), axis=0)]
    if X.shape[0] < 20:
        raise ValueError('too few samples for machine learning')
    if not os.path.isdir(args.output_dir):
        logger.info('create output directory: ' + args.output_dir)
        os.makedirs(args.output_dir)
    # read other input files
    logger.info('save class labels to: ' + os.path.join(args.output_dir, 'classes.txt'))
    pd.Series(y).to_csv(os.path.join(args.output_dir, 'classes.txt'), header=False, index=False)
    logger.info('save feature names to: ' + os.path.join(args.output_dir, 'feature_names.txt'))
    pd.Series(feature_names).to_csv(os.path.join(args.output_dir, 'feature_names.txt'), header=False, index=False)
    logger.info('save sample ids to: ' + os.path.join(args.output_dir, 'samples.txt'))
    pd.Series(sample_ids).to_csv(os.path.join(args.output_dir, 'samples.txt'), header=False, index=False)
    # read configuration file
    config = {}
    if args.config is not None:
        logger.info('read configuration file: ' + args.config)
        with open(args.config, 'r') as f:
            config = yaml.load(f)
    # overwride config with command-line pamameters
    argdict = vars(args)
    config.update(search_dict(argdict, (
        'zero_fraction_filter', 'zero_fraction_filter_params',
        'rpm_filter', 'rpm_filter_params',
        'rpkm_filter', 'rpkm_filter_params',
        'fold_change_filter', 'fold_change_filter_params',
        'diffexp_filter', 'diffexp_filter_params',
        'log_transform','log_transform_params',
        'scaler', 'scaler_params',
        'selector', 'selector_params', 'n_features_to_select',
        'classifier', 'classifier_params',
        'grid_search', 'grid_search_params'
    )))

    

    for key in ('rpkm_filter_params', 'rpm_filter_params', 'fold_change_filter_params',
        'zero_fraction_filter_params', 'log_transform_params', 'diffexp_filter_params',
        'scaler_params', 'classifier_params', 'selector_params', 'grid_search_params'):
        config[key] = parse_params(argdict[key])
    # set temp_dir for diffexp_selector
    if args.diffexp_filter:
        config['diffexp_filter_params']['temp_dir'] = os.path.join(args.output_dir, 'diffexp')
        logger.info('set temp_dir of diffexp_filter: {}'.format(config['diffexp_filter_params']['temp_dir']))
    logger.info('build combined estimator')
    estimator = CombinedEstimator(**config)

    logger.info('start cross-validation')
    collect_metrics = CollectMetrics()
    collect_predictions = CollectPredictions()
    collect_train_index = CollectTrainIndex()
    cv_callbacks = [collect_metrics, collect_predictions, collect_train_index]
    if args.selector is not None:
        feature_selection_matrix = FeatureSelectionMatrix()
        cv_callbacks.append(feature_selection_matrix)
    cv_params = parse_params(args.cv_params)
    if args.sample_weight is not None:
        if args.sample_weight == 'auto':
            sample_weight = 'auto'
        else:
            sample_weight = pd.read_table(args.sample_weight, header=None, index_col=0).iloc[:, 0]
    else:
        sample_weight = None
    _cross_validation(estimator, X, y, sample_weight=sample_weight, params=cv_params, callbacks=cv_callbacks)
    logger.info('collect_metrics:')
    #print(cv_callbacks[0].get_metrics())
    logger.info('fit estimator on full dataset')
    estimator.fit(X, y)
    logger.info('save final model to: ' + os.path.join(args.output_dir, 'final_model.pkl'))
    with open(os.path.join(args.output_dir, 'final_model.pkl'), 'wb') as f:
        pickle.dump(estimator, f)
    logger.info('classifier params: {}'.format(estimator.classifier_.get_params()))
    if args.selector is not None:
        feature_index = estimator.features_
        logger.info('number of selected features: {}'.format(feature_index.shape[0]))
        logger.info('save features to: ' + os.path.join(args.output_dir, 'features.txt'))
        pd.Series(feature_names[feature_index]).to_csv(os.path.join(args.output_dir, 'features.txt'), index=False, header=False)
        logger.info('save feature importances to: ' + os.path.join(args.output_dir, 'feature_importances.txt'))
        pd.Series(estimator.feature_importances_, index=feature_names[feature_index])\
            .to_csv(os.path.join(args.output_dir, 'feature_importances.txt'), sep='\t', header=False, index=True)
        #logger.info('save feature selection matrix to: ' + os.path.join(args.output_dir, 'feature_selection_matrix.txt'))
        #m = pd.DataFrame(feature_selection_matrix.get_matrix(), columns=feature_names)
        #m.columns.name = 'feature'
        #m.T.to_csv(os.path.join(args.output_dir, 'feature_selection_matrix.txt'), sep='\t', header=True, index=False)
    
    metrics = collect_metrics.get_metrics()
    for name in ('train', 'test'):
        logger.info('save metrics to: ' + os.path.join(args.output_dir, 'metrics.{}.txt'.format(name)))
        # if there are missing features, set metrics to NA
        metrics[name].to_csv(os.path.join(args.output_dir, 'metrics.{}.txt'.format(name)), header=True, index=True, na_rep='NA', sep='\t')

    logger.info('save cross-validation details to: ' + os.path.join(args.output_dir, 'cross_validation.h5'))
    with h5py.File(os.path.join(args.output_dir, 'cross_validation.h5'), 'w') as f:
        f.create_dataset('labels', data=y)
        f.create_dataset('predicted_labels', data=collect_predictions.get_pred_labels())
        f.create_dataset('predictions', data=collect_predictions.get_pred_probs())
        f.create_dataset('train_index', data=collect_train_index.get_train_index())
        if args.selector is not None:
            f.create_dataset('feature_selection', data=feature_selection_matrix.get_matrix())

@command_handler
def run_pipeline(args):
    from estimators2 import search_dict, CollectMetrics, CollectPredictions, CollectTrainIndex, FeatureSelectionMatrix,\
        CombinedEstimator, get_features_from_pipeline, parse_params
    from estimators2 import cross_validation as _cross_validation
    from sklearn.utils.class_weight import compute_sample_weight
    import pandas as pd
    import numpy as np
    import h5py
    import pickle
    import yaml

    logger.info('read configuration file: ' + args.config)
    with open(args.config, 'r') as f:
        config = yaml.load(f)
    # overwride config variables with command line arguments
    if args.matrix is not None:
        config['matrix'] = args.matrix
    if args.sample_classes is not None:
        config['sample_classes'] = args.sample_classes
    if args.positive_class is not None:
        config['positive_class'] = args.positive_class.split(',')
    if args.negative_class is not None:
        config['negative_class'] = args.negative_class.split(',')
    if args.features is not None:
        config['features'] = args.features

    # read input data matrix
    X, y, sample_ids, feature_names = read_data_matrix(
        config['matrix'], config['sample_classes'],
        **search_dict(config, ('features', 'transpose', 'positive_class', 'negative_class')))

    # fill missing features
    X = np.nan_to_num(X)
    #X = X[:, np.all(~np.isnan(X), axis=0)]
    if X.shape[0] < 20:
        raise ValueError('too few samples for machine learning')
    if not os.path.isdir(args.output_dir):
        logger.info('create output directory: ' + args.output_dir)
        os.makedirs(args.output_dir)
    # read other input files
    logger.info('save class labels to: ' + os.path.join(args.output_dir, 'classes.txt'))
    pd.Series(y).to_csv(os.path.join(args.output_dir, 'classes.txt'), header=False, index=False)
    logger.info('save feature names to: ' + os.path.join(args.output_dir, 'feature_names.txt'))
    pd.Series(feature_names).to_csv(os.path.join(args.output_dir, 'feature_names.txt'), header=False, index=False)
    logger.info('save sample ids to: ' + os.path.join(args.output_dir, 'samples.txt'))
    pd.Series(sample_ids).to_csv(os.path.join(args.output_dir, 'samples.txt'), header=False, index=False)

    # set temp_dir for diffexp_selector
    preprocess_steps = {}
    for step_dict in config['preprocess_steps']:
        step_tuple = tuple(step_dict.items())[0]
        preprocess_steps[step_tuple[0]] = step_tuple[1]
        preprocess_steps[step_tuple[0]]['params'] = step_tuple[1].get('params', {})
        if step_tuple[1]['name'] in ('DiffExpFilter', 'SIS'):
            temp_dir = os.path.join(args.output_dir, 'tmp')
            preprocess_steps[step_tuple[0]]['params']['temp_dir'] = temp_dir
            logger.info('set temp_dir of {} to {}'.format(step_tuple[1]['name'], temp_dir))
        
    logger.info('build combined estimator')
    estimator = CombinedEstimator(config)

    # add callbacks for cross-validation
    collect_metrics = CollectMetrics()
    collect_predictions = CollectPredictions()
    collect_train_index = CollectTrainIndex()
    cv_callbacks = [collect_metrics, collect_predictions, collect_train_index]
    # output feature selection if selector are found
    has_selector = False
    for step in preprocess_steps.values():
        if step['type'] == 'selector':
            has_selector = True
    if has_selector:
        logger.info('add cross-validation callback: FeatureSelectionMatrix')
        feature_selection_matrix = FeatureSelectionMatrix()
        cv_callbacks.append(feature_selection_matrix)
    # get feature weight
    if config.get('sample_weight') is not None:
        if config['sample_weight'] == 'balanced':
            logger.info('compute sample weight from class ratio')
            sample_weight = 'balanced'
        else:
            logger.info('read sample weight from file: ' + config['sample_weight'])
            sample_weight = pd.read_table(
                config['sample_weight'], header=None, index_col=0).iloc[:, 0]
    else:
        sample_weight = None
    logger.info('start cross-validation')
    _cross_validation(estimator, X, y, sample_weight=sample_weight, 
        params=config['cv_params'], callbacks=cv_callbacks)
    logger.info('collect_metrics:')
    #print(cv_callbacks[0].get_metrics())
    logger.info('fit estimator on full dataset')
    if config.get('sample_weight') == 'balanced':
        sample_weight = compute_sample_weight(class_weight='balanced', y=y)
    estimator.fit(X, y, sample_weight=sample_weight)
    logger.info('save final model to: ' + os.path.join(args.output_dir, 'final_model.pkl'))
    with open(os.path.join(args.output_dir, 'final_model.pkl'), 'wb') as f:
        pickle.dump(estimator, f)
    logger.info('classifier params: {}'.format(estimator.classifier_.get_params()))
    if has_selector:
        feature_index = estimator.features_
        logger.info('number of selected features: {}'.format(feature_index.shape[0]))
        logger.info('save features to: ' + os.path.join(args.output_dir, 'features.txt'))
        pd.Series(feature_names[feature_index]).to_csv(os.path.join(args.output_dir, 'features.txt'), index=False, header=False)
        logger.info('save feature importances to: ' + os.path.join(args.output_dir, 'feature_importances.txt'))
        pd.Series(estimator.feature_importances_, index=feature_names[feature_index])\
            .to_csv(os.path.join(args.output_dir, 'feature_importances.txt'), sep='\t', header=False, index=True)
        #logger.info('save feature selection matrix to: ' + os.path.join(args.output_dir, 'feature_selection_matrix.txt'))
        #m = pd.DataFrame(feature_selection_matrix.get_matrix(), columns=feature_names)
        #m.columns.name = 'feature'
        #m.T.to_csv(os.path.join(args.output_dir, 'feature_selection_matrix.txt'), sep='\t', header=True, index=False)
    
    metrics = collect_metrics.get_metrics()
    for name in ('train', 'test'):
        logger.info('save metrics to: ' + os.path.join(args.output_dir, 'metrics.{}.txt'.format(name)))
        # if there are missing features, set metrics to NA
        metrics[name].to_csv(os.path.join(args.output_dir, 'metrics.{}.txt'.format(name)), header=True, index=True, na_rep='NA', sep='\t')

    logger.info('save cross-validation details to: ' + os.path.join(args.output_dir, 'cross_validation.h5'))
    with h5py.File(os.path.join(args.output_dir, 'cross_validation.h5'), 'w') as f:
        f.create_dataset('labels', data=y)
        f.create_dataset('predicted_labels', data=collect_predictions.get_pred_labels())
        f.create_dataset('predictions', data=collect_predictions.get_pred_probs())
        f.create_dataset('train_index', data=collect_train_index.get_train_index())
        #print(feature_selection_matrix.get_matrix())
        if has_selector:
            f.create_dataset('feature_selection', data=feature_selection_matrix.get_matrix())

def predict(args):
    import pandas as pd
    import numpy as np

    logger.info('read data matrix: ' + matrix)
    X = pd.read_table(matrix, index_col=0, sep='\t')
    # transpose
    if transpose:
        logger.info('transpose feature matrix')
        X = X.T
    model_file = os.path.join(args.model_dir, 'final_model.pkl')
    logger.info('load model: ' + model_file)
    with open(model_file, 'rb') as f:
        model = pickle.load(f)
    logger.info('run model')
    predicted_scores = model.predict(X.values)
    logger.info('save predictions to file: ' + args.output_file)
    output_df = pd.DataFrame({'scores': predicted_scores[:, 1]})
    output_df.index = X.index.values
    output_df.index.name = 'sample_id'
    output_df.to_csv(args.output_file, sep='\t', na_rep='NA')

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
    g_input.add_argument('--negative-class', type=str, metavar='STRING',
        help='comma-separates list of sample classes to use as negative class')
    g_input.add_argument('--transpose', action='store_true', default=False,
        help='transpose the feature matrix')
    g_input.add_argument('--features', type=str, metavar='FILE',
        help='input file containing subset of feature names')
    g_input.add_argument('--config', '-c', type=str, metavar='FILE',
        help='configuration file of parameters in YAML format')

    g_filter = parser.add_argument_group('filter')
    g_filter.add_argument('--zero-fraction-filter', action='store_true')
    #g_filter.add_argument('--zero-fraction-filter-threshold', type=float, metavar='NUMBER')
    g_filter.add_argument('--zero-fraction-filter-params', type=str, metavar='STRING')
    g_filter.add_argument('--rpkm-filter', action='store_true')
    #g_filter.add_argument('--rpkm-filter-threshold', type=float, metavar='NUMBER')
    g_filter.add_argument('--rpkm-filter-params', type=str, metavar='STRING')
    g_filter.add_argument('--rpm-filter', action='store_true')
    #g_filter.add_argument('--rpm-filter-threshold', type=float, metavar='NUMBER')
    g_filter.add_argument('--rpm-filter-params', type=str, metavar='STRING')
    g_filter.add_argument('--fold-change-filter', action='store_true')
    #g_filter.add_argument('--fold-change-filter-direction', type=str, default='any', metavar='STRING')
    g_filter.add_argument('--fold-change-filter-params', type=str, metavar='STRING')
    g_filter.add_argument('--diffexp-filter', action='store_true')
    g_filter.add_argument('--diffexp-filter-params', type=str, metavar='STRING')

    g_scaler = parser.add_argument_group('scaler')
    g_scaler.add_argument('--log-transform', action='store_true')
    #g_scaler.add_argument('--log-transform-base', type=float, metavar='NUMBER')
    g_scaler.add_argument('--log-transform-params', type=str, metavar='STRING')
    g_scaler.add_argument('--scaler', type=str, metavar='NAME')
    g_scaler.add_argument('--scaler-params', type=str, metavar='STRING')

    g_select = parser.add_argument_group('feature_selection')
    g_select.add_argument('--selector', type=str, metavar='NAME')
    g_select.add_argument('--selector-params', type=str, metavar='STRING')
    g_select.add_argument('--n-features-to-select', type=int, metavar='INTEGER')

    g_classifier = parser.add_argument_group('classifier')
    g_classifier.add_argument('--classifier', type=str, metavar='NAME', default='random_forest')
    g_classifier.add_argument('--classifier-params', type=str, metavar='STRING')

    g_cv = parser.add_argument_group('cross_validation')
    g_cv.add_argument('--cv-params', type=str, metavar='STRING', nargs='?')
    g_cv.add_argument('--grid-search', action='store_true')
    g_cv.add_argument('--grid-search-params', type=str, metavar='STRING')

    g_misc= parser.add_argument_group('misc')
    g_misc.add_argument('--sample-weight', type=str, default='auto',
        help='''sample weight to balance classes. 
        Compute from data if set to "auto". 
        Can be a tab-separated file with two columns (no header): sample_id, weight.
        No sample weight if set to "none".''')
    
    g_output= parser.add_argument_group('output')
    g_output.add_argument('--output-dir', '-o', type=str, metavar='DIR', 
        required=True, help='output directory')
    
    parser = subparsers.add_parser('run_pipeline')
    parser.add_argument('--matrix', '-i', type=str, metavar='FILE', required=True,
        help='input feature matrix (rows are samples and columns are features')
    parser.add_argument('--sample-classes', type=str, metavar='FILE', required=True,
        help='input file containing sample classes with 2 columns: sample_id, sample_class')
    parser.add_argument('--positive-class', type=str, metavar='STRING',
        help='comma-separated list of sample classes to use as positive class')
    parser.add_argument('--negative-class', type=str, metavar='STRING',
        help='comma-separates list of sample classes to use as negative class')
    parser.add_argument('--features', type=str, metavar='FILE',
        help='input file containing subset of feature names')
    parser.add_argument('--config', '-c', type=str, metavar='FILE', required=True,
        help='configuration file of parameters in YAML format')
    parser.add_argument('--output-dir', '-o', type=str, metavar='DIR', 
        required=True, help='output directory')
    parser.add_argument('--log-level', type=str, default='INFO',
        help='logging level')

    parser = subparsers.add_parser('predict')
    parser.add_argument('--matrix', '-i', type=str, metavar='FILE', required=True,
        help='input feature matrix (rows are samples and columns are features')
    parser.add_argument('--model-dir', '-m', type=str, metavar='DIR', required=True,
        help='directory generated by feature_selection')
    parser.add_argument('--transpose', action='store_true', default=False,
        help='transpose the input matrix')
    parser.add_argument('--output-file', '-o', type=str, metavar='DIR', required=True,
        help='output file')
    
    parser = subparsers.add_parser('generate_cv_splits')
    parser.add_argument('--samples', '-i', type=str, required=True,
        help='input file to determine number of samples by n_lines - 1')
    parser.add_argument('--splitter', type=str, default='KFold')
    parser.add_argument('--spliter-params', type=str,
        help='parameters of the splitter')
    parser.add_argument('--output-file', '-o', type=str, required=True,
        help='output file')

    args = main_parser.parse_args()
    if args.command is None:
        print('Errror: missing command', file=sys.stdout)
        main_parser.print_help()
        sys.exit(1)

    import pandas as pd
    import numpy as np

    command_handlers.get(args.command)(args)