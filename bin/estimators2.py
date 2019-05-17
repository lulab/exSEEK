import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin, MetaEstimatorMixin, is_classifier
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC, SVC
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve, \
    precision_score, recall_score, f1_score, \
    average_precision_score
from sklearn.feature_selection.base import SelectorMixin
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler, MaxAbsScaler
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.model_selection import GridSearchCV, check_cv
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_X_y
from abc import ABC, ABCMeta, abstractmethod
from tqdm import tqdm
import pickle
import json
import os
import subprocess
import shutil
import inspect
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm
import logging
from copy import deepcopy
from utils import search_dict, get_feature_importances, function_has_arg
from feature_selectors import get_selector

logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger(__name__)

def parse_params(s):
    '''Parse param dict from string
    
    Returns:
        params: dict
        If s is None, return empty dict
        If s is in JSON format, return parsed object
    '''
    if s is not None:
        params = s
        try:
            params = json.loads(s)
        except json.JSONDecodeError:
            pass
        finally:
            return params
    else:
        return {}



def predict_proba(estimator, X):
    try:
        proba = estimator.predict_proba(X)
    except AttributeError:
        proba = estimator.decision_function(X)
    return proba

def get_scorer(scoring):
    '''Get scoring function from string

    Parameters:
        scoring: str
            choices: roc_auc, accuracy
    
    Returns:
        score_func: function(y_true, y_pred)
    '''
    if scoring == 'roc_auc':
        return roc_auc_score
    elif scoring == 'accuracy':
        return accuracy_score
    else:
        raise ValueError('unknonwn scoring: {}'.format(scoring))

def classification_scores(y_true, y_pred_labels, y_pred_probs):
    scores = {
        'roc_auc': roc_auc_score(y_true, y_pred_probs),
        'average_precision': average_precision_score(y_true, y_pred_probs),
        'accuracy': accuracy_score(y_true, y_pred_labels),
        'precision': precision_score(y_true, y_pred_labels),
        'recall': recall_score(y_true, y_pred_labels),
        'f1_score': f1_score(y_true, y_pred_labels)
    }
    return scores

def get_classifier(name, **params):
    '''Get scoring function from string

    Parameters:
        name: str
            name of the clasifier

        params: keyword arguments
            extra parameters for the classifier
    
    Returns:
        estimator: object
            a BaseEstimator object
    '''
    if name == 'LogisticRegression':
        return LogisticRegression(**search_dict(params, 
            ('penalty', 'dual', 'C', 'tol', 'fit_intercept', 'solver',
                'class_weight', 'max_iter', 'n_jobs', 'random_state', 'verbose')))
    elif name == 'LogisticRegressionL1':
        return LogisticRegression(penalty='l1', **search_dict(params, 
            ('dual', 'C', 'tol', 'fit_intercept', 'solver',
            'class_weight', 'max_iter', 'n_jobs', 'random_state', 'verbose')))
    elif name == 'LogisticRegressionL2':
        return LogisticRegression(penalty='l2', **search_dict(params, 
            ('dual', 'C', 'tol', 'fit_intercept', 'solver',
            'class_weight', 'max_iter', 'n_jobs', 'random_state', 'verbose')))
    elif name == 'RandomForestClassifier':
        return RandomForestClassifier(**search_dict(params,
            ('n_estimators', 'criterion', 'max_depth', 'min_samples_split', 'min_samples_leaf', 
             'min_weight_fraction_leaf', 'max_features', 'max_leaf_nodes', 
             'min_impurity_decrease', 'min_impurity_split', 'oob_score',
             'n_jobs', 'verbose', 'random_state', 'class_weight')))
    elif name == 'LinearSVC':
        return LinearSVC(**search_dict(params,
            ('penalty', 'loss', 'dual', 'tol', 'C', 'fit_intercept', 
             'intercept_scaling', 'class_weight', 'verbose',
             'random_state', 'max_iter')))
    elif name == 'SVC':
        return SVC(**search_dict(params,
            ('penalty', 'loss', 'dual', 'tol', 'C', 'fit_intercept', 'gamma',
             'intercept_scaling', 'class_weight', 'verbose',
             'random_state', 'max_iter')))
    elif name == 'DecisionTreeClassifier':
        return DecisionTreeClassifier(**search_dict(params,
            ('criterion', 'splitter', 'max_depth', 'min_samples_split', 'min_samples_leaf',
             'min_weight_fraction_leaf', 'max_features', 'max_leaf_nodes', 'min_impurity_decrease',
             'min_impurity_split')))
    elif name == 'ExtraTreesClassifier':
        return ExtraTreesClassifier(**search_dict(params,
            ('n_estimators', 'criterion', 'max_depth', 'min_samples_split', 'min_samples_leaf', 
             'min_weight_fraction_leaf', 'max_features', 'max_leaf_nodes', 
             'min_impurity_decrease', 'min_impurity_split', 'oob_score',
             'n_jobs', 'verbose', 'random_state', 'class_weight')))
    elif name == 'MLPClassifier':
        from sklearn.neural_network import MLPClassifier
        return MLPClassifier(**search_dict(params,
            ('hidden_layer_sizes', 'activation', 'solver',
             'alpha', 'batch_size', 'learning_rate', 'max_iter')))
    elif name == 'SGDClassifier':
        from sklearn.linear_model import SGDClassifier
        return SGDClassifier(**search_dict(params,
            ('loss', 'penalty', 'alpha', 'l1_ratio', 'fit_intercept',
            'max_iter', 'tol', 'epsilon')))
    else:
        raise ValueError('unknown classifier: {}'.format(name))

def get_transformer(name, **params):
    if name == 'LogTransform':
        return LogTransform(**params)
    else:
        raise ValueError('unknown transformer: {}'.format(name))

class FileSplitter(object):
    '''Splitter that read train and test index from given file
    Input file is a tab-separated file. Each row is a train-test split.
    The number of columns equals to the number of samples. 
    1's indicate training samples and 0's indicate test samples
    
    '''
    def __init__(self, filename):
        self.filename = filename
        self.cv_matrix = np.loadtxt(filename, dtype='int', delimiter='\t')

    def split(self, X, y=None):
        for i in range(self.cv_matrix.shape[0]):
            train_index = np.nonzero(self.cv_matrix[i])[0]
            test_index = np.nonzero(~self.cv_matrix[i])[0]
            yield train_index, test_index

def get_splitter(random_state=None, **params):
    '''Get cross-validation index generator

    Parameters:
        random_state: int or RandomState object
            seed for random number generator 
        
        name: str
            name of the splitter

        params: keyword arguments
            extra parameters for the classifier
    
    Returns:
        estimator: object
            a BaseEstimator object
    '''
    from sklearn.model_selection import KFold, StratifiedKFold, ShuffleSplit, LeaveOneOut, \
        RepeatedKFold, RepeatedStratifiedKFold, LeaveOneOut, StratifiedShuffleSplit

    splitter = params.get('splitter')
    if splitter is None:
        return check_cv(**params)
    if splitter == 'KFold':
        from sklearn.model_selection import KFold
        return KFold(random_state=random_state, **search_dict(params, ('n_splits', 'shuffle')))
    elif splitter == 'StratifiedKFold':
        from sklearn.model_selection import StratifiedKFold
        return StratifiedKFold(random_state=random_state, **search_dict(params, ('n_splits', 'shuffle')))
    elif splitter == 'RepeatedStratifiedKFold':
        from sklearn.model_selection import RepeatedStratifiedKFold
        return RepeatedStratifiedKFold(random_state=random_state, **search_dict(params, ('n_splits', 'n_repeats')))
    elif splitter == 'ShuffleSplit':
        from sklearn.model_selection import ShuffleSplit
        return ShuffleSplit(random_state=random_state, **search_dict(params, ('n_splits', 'test_size', 'train_size')))
    elif splitter == 'StratifiedShuffleSplit':
        from sklearn.model_selection import StratifiedShuffleSplit
        return StratifiedShuffleSplit(random_state=random_state, **search_dict(params, ('n_splits', 'test_size', 'train_size')))
    elif splitter == 'LeaveOneOut':
        from sklearn.model_selection import LeaveOneOut
        return LeaveOneOut()
    elif splitter == 'FileSplitter':
        return UserFileSplitter(**search_dict(params, 'filename'))
    else:
        raise ValueError('unknown splitter: {}'.format(splitter))
        

def get_score_function(estimator):
    '''Get method of an estimator that predict a continous score for each sample
    '''
    if hasattr(estimator, 'predict_proba'):
        return estimator.predict_proba
    elif hasattr(estimator, 'decision_function'):
        return estimator.decision_function
    else:
        raise ValueError('the estimator should either have decision_function() method or predict_proba() method')

def get_scaler(name, **params):
    if name == 'StandardScaler':
        return StandardScaler(**search_dict(params, ('with_mean', 'with_std', 'copy')))
    elif name == 'RobustScaler':
        return RobustScaler(**search_dict(params, ('with_centering', 'with_scaling', 'quantile_range', 'copy')))
    elif name == 'MinMaxScaler':
        return MinMaxScaler(**search_dict(params, ('feature_range', 'copy')))
    elif name == 'MaxAbs':
        return MaxAbsScaler(**search_dict(params, ('copy',)))
    elif name == 'LogTransform':
        return LogTransform(**search_dict(params, ('base', 'pseudo_count')))


class LogTransform(BaseEstimator, TransformerMixin):
    '''Transform features by applying logarithm function

    Parameters:
    ----------

    base: float
        The logarithm base. Natural logarithm is used if set to None.
    
    pseudo_count: float
        Pseudo-count added to the original matrix.
    '''
    def __init__(self, base=None, pseudo_count=0.01):
        self.base = None
        self.pseudo_count = pseudo_count
    
    def fit(self, X, y=None, **kwargs):
        return self
    
    def transform(self, X, y=None):
        if self.pseudo_count != 0:
            X = X + self.pseudo_count
        if self.base is None:
            return np.log(X)
        elif self.base == 2:
            return np.log2(X)
        elif self.base == 10:
            return np.log10(X)
        else:
            return np.log(X)/np.log(self.base)
    
    def inverse_transform(self, X, y=None):
        if self.base is None:
            X = np.exp(X)
        else:
            X = np.power(X, self.base)
        if self.pseudo_count != 0:
            X -= self.pseudo_count
        return X

def get_features_from_pipeline(pipeline, n_features):
    X = np.arange(n_features).reshape((1, -1))
    for name, step in pipeline.named_steps.items():
        if isinstance(step, SelectorMixin):
            X = step.transform(X)
    return np.ravel(X)

class CVCallback(ABC):
    @abstractmethod
    def __call__(self, estimator, X, y, y_pred_labels, y_pred_probs, train_index, test_index):
        pass

class CollectMetrics(CVCallback):
    def __init__(self, scoring='roc_auc', classifier='classifier', has_missing_features=False):
        self.scoring = scoring
        self.metrics = {'train': [], 'test': []}
        self.classifier = classifier
        self.has_missing_features = has_missing_features
    
    def __call__(self, estimator, X, y, y_pred_labels, y_pred_probs, train_index, test_index):
        self.metrics['train'].append(classification_scores(
            y[train_index], y_pred_labels[train_index], y_pred_probs[train_index]))
        self.metrics['test'].append(classification_scores(
            y[test_index], y_pred_labels[test_index], y_pred_probs[test_index]))
    
    def get_metrics(self):
        for name in ('train', 'test'):
            if isinstance(self.metrics[name], list):
                self.metrics[name] = pd.DataFrame.from_records(self.metrics[name])
                self.metrics[name].index.name = 'split'
                if self.has_missing_features:
                    self.metrics[name][:] = np.nan
        return self.metrics

class CollectPredictions(CVCallback):
    def __init__(self):
        self.pred_labels = []
        self.pred_probs = []
    
    def __call__(self, estimator, X, y, y_pred_labels, y_pred_probs, train_index, test_index):
        self.pred_labels.append(np.ravel(y_pred_labels))
        self.pred_probs.append(y_pred_probs)
    
    def get_pred_labels(self):
        if isinstance(self.pred_labels, list):
            self.pred_labels = np.vstack(self.pred_labels)
        return self.pred_labels
    
    def get_pred_probs(self):
        if isinstance(self.pred_probs, list):
            self.pred_probs = np.vstack(self.pred_probs)
        return self.pred_probs

class FeatureSelectionMatrix(CVCallback):
    def __init__(self, selector='selector'):
        self.matrix = []
        self.selector = selector
    
    def __call__(self, estimator, X, y, y_pred_labels, y_pred_probs, train_index, test_index):
        support = np.zeros(X.shape[1], dtype='bool')
        support[estimator.features_] = True
        self.matrix.append(support)
    
    def get_matrix(self):
        if isinstance(self.matrix, list):
            self.matrix = np.vstack(self.matrix)
        return self.matrix

class CollectTrainIndex(CVCallback):
    def __init__(self):
        self.train_index = []
    
    def __call__(self, estimator, X, y, y_pred_labels, y_pred_probs, train_index, test_index):
        ind = np.zeros(X.shape[0], dtype='bool')
        ind[train_index] = True
        self.train_index.append(ind)
    
    def get_train_index(self):
        if isinstance(self.train_index, list):
            self.train_index = np.vstack(self.train_index)
        return self.train_index

class CombinedEstimator1(BaseEstimator, MetaEstimatorMixin):
    def __init__(self,
        zero_fraction_filter=False,
        zero_fraction_filter_params=None,
        fold_change_filter=False,
        fold_change_filter_params=None,
        rpkm_filter=False,
        rpkm_filter_params=None,
        rpm_filter=False,
        rpm_filter_params=None,
        log_transform=False,
        log_transform_params=None,
        diffexp_filter=None,
        diffexp_filter_params=None,
        scaler=None,
        scaler_params=None,
        selector=None,
        selector_params=None,
        n_features_to_select=None,
        classifier=None,
        classifier_params=None,
        grid_search=False,
        grid_search_params=None,
        **kwargs):
        self.zero_fraction_filter = zero_fraction_filter
        self.zero_fraction_filter_params = zero_fraction_filter_params
        
        self.rpkm_filter = rpkm_filter
        self.rpkm_filter_params = rpkm_filter_params
        
        self.rpm_filter = rpm_filter
        self.rpm_filter_params = rpm_filter_params
        
        self.fold_change_filter = fold_change_filter
        self.fold_change_filter_params = fold_change_filter_params
        
        self.log_transform = log_transform
        self.log_transform_params = log_transform_params
        self.diffexp_filter = diffexp_filter
        self.diffexp_filter_params = diffexp_filter_params
        self.scaler = scaler
        self.scaler_params = scaler_params
        self.grid_search = grid_search
        self.grid_search_params = grid_search_params
        self.classifier = classifier
        self.classifier_params = classifier_params
        self.selector = selector
        self.selector_params = selector_params
        self.n_features_to_select = n_features_to_select
        
    """
    def __init__(self, config):
        self.config = config
    """
    @staticmethod
    def get_gene_lengths_from_feature_names(feature_names):
        feature_info = pd.Series(feature_names).str.split('|', expand=True)
        feature_info.columns = ['gene_id', 'gene_type', 'gene_name', 'feature_id', 'transcript_id', 'start', 'end']
        feature_info['start'] = feature_info['start'].astype('int')
        feature_info['end'] = feature_info['end'].astype('int')
        feature_info['length'] = feature_info['end'] - feature_info['start']
        return feature_info['length'].values
    
    def fit(self, X, y=None, sample_weight=None):
        self.preprocess_steps_ = []
        if self.zero_fraction_filter:
            logger.debug('add zero_fraction_filter with parameters: {}'.format(self.zero_fraction_filter_params))
            self.preprocess_steps_.append(('zero_fraction_filter', 
                get_selector('zero_fraction_filter', **self.zero_fraction_filter_params)))
        '''
        if self.rpkm_filter:
            logger.debug('add rpkm_filter with parameters: {}'.format(self.rpkm_filter_params))
            if self.feature_names is None:
                raise ValueError('feature_names is required for rpkm_filter')
            gene_lengths = self.get_gene_lengths_from_feature_names(feature_names)
            step = get_selector('rpkm_filter', **rpkm_filter_params)
            step.set_gene_lengths(gene_lengths)
            preprocess_steps.append(('rpkm_filter', step))
        '''
        if self.rpm_filter:
            logger.debug('add rpm_filter with parameters: {}'.format(self.rpm_filter_params))
            self.preprocess_steps_.append(('rpm_filter', get_selector('rpm_filter', **self.rpkm_filter_params)))
        if self.fold_change_filter:
            logger.debug('add fold_change_filter with parameters: {}'.format(self.fold_change_filter_params))
            self.preprocess_steps_.append(('fold_change_filter',
                get_selector('fold_change_filter', **self.fold_change_filter_params)))
        if self.diffexp_filter:
            logger.debug('add diffexp_filter with parameters: {}'.format(self.diffexp_filter_params))
            self.preprocess_steps_.append(('diffexp_filter',
                get_selector('diffexp_filter', **self.diffexp_filter_params)))
        if self.log_transform:
            logger.debug('add log_transform with parameters: {}'.format(self.log_transform_params))
            self.preprocess_steps_.append(('log_transform', 
                get_scaler('log_transform', **self.log_transform_params)))
        if self.scaler is not None:
            logger.debug('add scaler "{}" with parameters: {}'.format(self.scaler, self.scaler_params))
            self.preprocess_steps_.append(('scaler', 
                get_scaler(self.scaler, **self.scaler_params)))
        # preprocess features
        X_new = X
        self.features_ = np.arange(X.shape[1])
        for name, step in self.preprocess_steps_:
            X_new = step.fit_transform(X_new, y)
            setattr(self, name + '_', step)
            if isinstance(step, SelectorMixin):
                self.features_ = self.features_[step.get_support()]
        logger.debug('add classifier "{}" with parameters: {}'.format(self.classifier, self.classifier_params))
        self.classifier_ = get_classifier(self.classifier, **self.classifier_params)
        # grid search for hyper-parameters
        if self.grid_search:
            logger.debug('add grid_search with parameters: {}'.format(self.grid_search_params))
            grid_search_params = deepcopy(self.grid_search_params)
            if 'cv' in grid_search_params:
                grid_search_params['cv'] = get_splitter(**grid_search_params['cv'])
            grid_search_params['param_grid'] = grid_search_params['param_grid'][self.classifier]
            self.grid_search_ = GridSearchCV(estimator=self.classifier_,
                **search_dict(grid_search_params, ('param_grid', 'scoring', 'cv', 
                'fit_params', 'verbose', 'return_train_score', 'error_score', 'iid')))
            self.grid_search_.fit(X_new, y, sample_weight=sample_weight)
            self.classfier_ = self.grid_search_.best_estimator_
            self.best_classifier_params_ = self.grid_search_.best_params_
            self.classifier_.set_params(**self.grid_search_.best_params_)
            #logger.info('best params: {}'.format(self.grid_search_.best_params_))
            #logger.info('mean test score: {}'.format(self.grid_search_.cv_results_['mean_test_score']))
        # feature selection
        if self.selector:
            logger.debug('add selector "{}" with parameters: {}'.format(self.selector, self.selector_params))
            logger.debug('number of features to select: {}'.format(self.n_features_to_select))
            # classifier for feature selection wrapper
            selector_classifier = None
            if 'classifier' in self.selector_params:
                selector_classifier = get_classifier(self.selector_params['classifier'], 
                    **self.selector_params['classifier_params'])
            else:
                selector_classifier = self.classifier_
            self.selector_ = get_selector(self.selector, estimator=selector_classifier, 
                n_features_to_select=self.n_features_to_select, **self.selector_params)
            X_new = self.selector_.fit_transform(X_new, y)
            self.features_ = self.features_[self.selector_.get_support()]
        # refit the classifier with selected features
        self.classifier_.fit(X_new, y, sample_weight=sample_weight)
        # set feature importances
        self.feature_importances_ = get_feature_importances(self.classifier_)
        return self
    
    def transform(self, X, y=None):
        check_is_fitted(self, 'classifier_')
        for name, step in self.preprocess_steps_:
            X = step.transform(X)
        if self.selector is not None:
            X = self.selector_.transform(X)
        return X
    
    def predict(self, X):
        X = self.transform(X)
        return self.classifier_.predict(X)
    
    def predict_proba(self, X):
        X = self.transform(X)
        try:
            proba = self.classifier_.predict_proba(X)[:, 1]
        except AttributeError:
            proba = self.classifier_.decision_function(X)
        return proba

def cross_validation(estimator, X, y, sample_weight='balanced', params=None, callbacks=None):
    splitter = get_splitter(**params)
    logger.debug('start cross-validation')
    logger.debug('cross-validation parameters: {}'.format(params))
    logger.debug('number of cross-validation splits: {}'.format(splitter.get_n_splits(X, y)))
    pbar = tqdm(unit='split', total=splitter.get_n_splits(X, y))
    for index in splitter.split(X, y):
        train_index, test_index = index
        estimator = clone(estimator)
        sample_weight_ = sample_weight
        if sample_weight == 'balanced':
            sample_weight_ = compute_sample_weight(class_weight='balanced', y=y[train_index])
        else:
            sample_weight_ = sample_weight[train_index]
        estimator.fit(X[train_index], y[train_index], sample_weight=sample_weight_)
        y_pred_labels = estimator.predict(X)
        y_pred_probs = estimator.predict_proba(X)
        for callback in callbacks:
            callback(estimator, X, y, y_pred_labels, y_pred_probs, train_index, test_index)
        pbar.update(1)
    pbar.close()

class CombinedEstimator(BaseEstimator, MetaEstimatorMixin):
    def __init__(self, config):
        self.config = config
    
    def fit(self, X, y=None, sample_weight=None):
        # use all features in the initial step
        self.features_ = np.arange(X.shape[1]).reshape((1, -1))
        X_new = X
        # preprocess steps
        self.preprocess_steps = []
        for step_dict in self.config['preprocess_steps']:
            step_name, step = list(step_dict.items())[0]
            if not step.get('enabled', False):
                continue
            if step['type'] == 'scaler':
                logger.debug('add scaler: {}.{}'.format(step_name, step['name']))
                preprocessor = get_scaler(step['name'], **step['params'])
            elif step['type'] == 'selector':
                # build classifier for wrapper-based selector
                selector_params = deepcopy(step['params'])
                # wrapper-based selector needs a classifier
                logger.debug('add selector: {}.{}'.format(step_name, step['name']))
                if selector_params.get('classifier') is not None:
                    logger.debug('get internal classifier: {}'.format(selector_params['classifier']))
                    selector_params['classifier_params'] = selector_params.get('classifier_params', {})
                    classifier = get_classifier(selector_params['classifier'], 
                        **selector_params['classifier_params'])
                    del selector_params['classifier']
                    del selector_params['classifier_params']
                    # optimize hyper-parameters by grid search
                    #logger.debug(selector_params)
                    if selector_params.get('grid_search', False):
                        logger.debug('grid search for internal classifier')
                        grid_search_params = deepcopy(selector_params['grid_search_params'])
                        grid_search_params['cv'] = get_splitter(
                            **grid_search_params['cv'])
                        grid_search = GridSearchCV(classifier, **grid_search_params)
                        grid_search.fit(X, y, sample_weight=sample_weight)
                        logger.debug('optimized hyper-parameters for internal classifier: {}'.format(grid_search.best_params_))
                        classifier = grid_search.best_estimator_
                        classifier.set_params(**grid_search.best_params_)
                        del grid_search
                        del selector_params['grid_search']
                        del selector_params['grid_search_params']
                elif step['type'] == 'transformer':
                    logger.debug('add transformer: {}.{}'.format(step_name, step['name']))
                    preprocessor = get_transformer(step['name'], **step['params'])
                else:
                    classifier = None
                preprocessor = get_selector(step['name'], classifier, 
                    n_features_to_select=self.config.get('n_features_to_select'), **selector_params)
            elif step['type'] == 'transformer':
                preprocessor = get_transformer(step['name'], **step['params'])
            else:
                raise ValueError('invalid preprocess step type: {}'.format(step['type']))
            # run the preprocessor
            if function_has_arg(preprocessor.fit, 'sample_weight'):
                preprocessor.fit(X_new, y, sample_weight=sample_weight)
            else:
                preprocessor.fit(X_new, y)
            X_new = preprocessor.transform(X_new)
            if step['type'] == 'selector':
                self.features_ = preprocessor.transform(self.features_)
                #self.features_ = self.features_[preprocessor.get_support()]
            # save the preprocessor
            self.preprocess_steps.append(preprocessor)
        # flatten feature indices
        self.features_ = self.features_.flatten()
        logger.debug('number of selected features: {}'.format(self.features_.shape[0]))
        # build classifier
        logger.debug('add classifier {}'.format(self.config['classifier']))
        classifier_params = self.config.get('classifier_params', {})
        self.classifier_ = get_classifier(self.config['classifier'], **classifier_params)
        # grid search for hyper-parameters
        if self.config.get('grid_search', False):
            logger.debug('grid search for classifier')
            grid_search_params = deepcopy(self.config['grid_search_params'])
            # get cross-validation splitter
            if grid_search_params.get('cv') is not None:
                grid_search_params['cv'] = get_splitter(**grid_search_params['cv'])
            grid_search_params['param_grid'] = grid_search_params['param_grid']
            self.grid_search_ = GridSearchCV(estimator=self.classifier_,
                **grid_search_params)
            if function_has_arg(self.classifier_.fit, 'sample_weight'):
                self.grid_search_.fit(X_new, y, sample_weight=sample_weight)
            else:
                self.grid_search_.fit(X_new, y)
            self.classifier_ = self.grid_search_.best_estimator_
            self.best_classifier_params_ = self.grid_search_.best_params_
            self.classifier_.set_params(**self.grid_search_.best_params_)
            logger.debug('best params: {}'.format(self.grid_search_.best_params_))
            #logger.info('mean test score: {}'.format(self.grid_search_.cv_results_['mean_test_score']))

        # refit the classifier with selected features
        if function_has_arg(self.classifier_.fit, 'sample_weight'):
            self.classifier_.fit(X_new, y, sample_weight=sample_weight)
        else:
            self.classifier_.fit(X_new, y)
        # set feature importances
        self.feature_importances_ = get_feature_importances(self.classifier_)
        return self
    
    def transform(self, X, y=None):
        check_is_fitted(self, 'classifier_')
        for step in self.preprocess_steps:
            X = step.transform(X)
        return X
    
    def predict(self, X):
        X = self.transform(X)
        return self.classifier_.predict(X)
    
    def predict_proba(self, X):
        X = self.transform(X)
        try:
            proba = self.classifier_.predict_proba(X)[:, 1]
        except AttributeError:
            proba = self.classifier_.decision_function(X)
        return proba