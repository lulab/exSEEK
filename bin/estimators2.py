import numpy as np

from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin, MetaEstimatorMixin, is_classifier
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import LinearSVC
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve
from sklearn.feature_selection.base import SelectorMixin
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler, MaxAbsScaler
from sklearn.model_selection import KFold, StratifiedKFold, ShuffleSplit, LeaveOneOut, \
        RepeatedKFold, RepeatedStratifiedKFold, LeaveOneOut, StratifiedShuffleSplit
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.model_selection import GridSearchCV, check_cv
from sklearn.feature_selection import RFE, RFECV
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_X_y
from abc import ABC, abstractmethod
import pandas as pd

from tqdm import tqdm
import pickle
import json
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

def parse_params(s):
    if s is not None:
        params = s
        try:
            params = json.loads(s)
        except json.JSONDecodeError as e:
            pass
        finally:
            return params
    else:
        return {}

def search_dict(data, keys):
    '''Intersect given keys with a dict and return the subset
    '''
    return {key:data[key] for key in keys if key in data}

def get_scorer(scoring):
    if scoring == 'roc_auc':
        return roc_auc_score
    elif scoring == 'accuracy':
        return accuracy_score
    else:
        raise ValueError('unknonwn scoring: {}'.format(scoring))

def get_classifier(name, **params):
    if name == 'logistic_regression':
        return LogisticRegression(**search_dict(params, 
            ('penalty', 'dual', 'C', 'tol', 'fit_intercept', 
                'class_weight', 'max_iter', 'n_jobs', 'random_state', 'verbose')))
    elif name == 'random_forest':
        return RandomForestClassifier(**search_dict(params,
            ('n_estimators', 'criterion', 'max_depth', 'min_samples_split', 'min_samples_leaf', 
             'min_weight_fraction_leaf', 'max_fetures', 'max_leaf_nodes', 
             'min_impurity_decrease', 'min_impurity_split', 'oob_score',
             'n_jobs', 'verbose', 'random_state', 'class_weight')))
    elif name == 'linear_svm':
        return LinearSVC(**search_dict(params,
            ('penalty', 'loss', 'dual', 'tol', 'C', 'fit_intercept', 
             'intercept_scaling', 'class_weight', 'verbose',
             'random_state', 'max_iter')))
    else:
        raise ValueError('unknown classifier: {}'.format(name))

def get_selector(name, estimator=None, **params):
    if name == 'robust':
        return RobustSelector(estimator, **search_dict(params,
         ('cv', 'n_features_to_select', 'verbose')))
    elif name == 'rfe':
        return RFE(estimator, **search_dict(params, 
        ('n_features_to_select', 'step', 'verbose')))
    elif name == 'rfecv':
        return RFECV(estimator, **search_dict(params,
         ('n_features_to_select', 'step', 'cv', 'verbose')))
    elif name == 'fold_change':
        return FoldChangeFilter(**search_dict(params,
        ('threshold', 'direction', 'below', 'pseudo_count')))
    elif name == 'zero_fraction':
        return ZeroFractionFilter(**search_dict(params,
        ('max_zero_fractino',)))
    elif name == 'rpkm_filter':
        return RpkmFilter(**search_dict(params,
        ('threshold',)))
    else:
        raise ValueError('unknown selector: {}'.format(name))

def get_splitter(**params):
    splitter = params.get('splitter')
    if splitter is None:
        return check_cv(**params)
    if splitter == 'kfold':
        return KFold(**search_dict(params, ('n_splits', 'shuffle', 'random_state')))
    elif splitter == 'stratified_kfold':
        return StratifiedKFold(**search_dict('n_splits', 'shuffle', 'random_state'))
    elif splitter == 'repeated_stratified_kfold':
        return RepeatedStratifiedKFold(**search_dict(params, ('n_splits', 'n_repeats', 'random_state')))
    elif splitter == 'shuffle_split':
        return ShuffleSplit(**search_dict(params, ('n_splits', 'test_size', 'train_size', 'random_state')))
    elif splitter == 'stratified_shuffle_split':
        return StratifiedShuffleSplit(**search_dict(params, ('n_splits', 'test_size', 'train_size', 'random_state')))
    elif splitter == 'leave_one_out':
        return LeaveOneOut()
    else:
        raise ValueError('unknown splitter: {}'.format(splitter))
        
def get_feature_importances(estimator):
    '''Get feature importance attribute of an estimator
    '''
    if hasattr(estimator, 'coef_'):
        return np.ravel(np.abs(estimator.coef_))
    elif hasattr(estimator, 'feature_importances_'):
        return np.ravel(estimator.feature_importances_)
    elif isinstance(estimator, RFE):
        ranking = estimator.ranking_.astype('float') - 1
        return np.ravel(1.0 - ranking/ranking.max())
    else:
        raise ValueError('the estimator should have either coef_ or feature_importances_ attribute')

def get_feature_ranking(feature_importances):
    '''Calculate ranking from feature importances
    Feature importances are sorted in descendig order and then converted to ranks
    Smaller values indicate higher importance

    Parameters
    ----------
    arrays: list of array-like objects
        feature importances
    
    Returns
    -------
    arrays: array-like
        Feature ranking
    '''
    ranking = np.zeros(len(feature_importances), dtype='int')
    ranking[np.argsort(-feature_importances)] = np.arange(len(feature_importances))
    return ranking

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
    if name == 'zscore':
        return StandardScaler(**search_dict(params, ('with_mean', 'with_std', 'copy')))
    elif name == 'robust':
        return RobustScaler(**search_dict(params, ('with_centering', 'with_scaling', 'quantile_range', 'copy')))
    elif name == 'min_max':
        return MinMaxScaler(**search_dict(params, ('feature_range', 'copy')))
    elif name == 'max_abs':
        return MaxAbsScaler(**search_dict(params, ('copy',)))
    elif name == 'log_transform':
        return LogTransform(**search_dict(params, ('base', 'pseudo_count')))

def bootstrap(*arrays):
    '''Perform bootstrap resampling

    Parameters
    ----------
    arrays: list of array-like objects
        Data to resample. Bootstrap resampling is performed on the first dimension of each array
    
    Returns
    -------
    arrays: tuple of array-like objects
        Resampled data
    '''
    n_samples = arrays[0].shape[0]
    if not all(n_samples == a.shape[0] for a in arrays):
        raise ValueError('the first dimension of all arrays should be equal')
    indices = np.random.randint(n_samples, size=n_samples)
    return (np.take(a, indices, axis=0) for a in arrays)

def jackknife(*arrays, remove=1, indices=None):
    '''Perform jackknife resampling

    Parameters
    ----------
    arrays: list of array-like objects
        Data to resample. Bootstrap resampling is performed on the first dimension of each array
    
    remove: int or float
        If `remove` is a integer, remove that number of samples.
        If `remove` is a float, remove that fraction of samples.

    indices: array-like
        Indices of samples to remove. The `remove` parameter will be ignored.

    Returns
    -------
    arrays: tuple of array-like objects
        Resampled data
    '''
    n_samples = arrays[0].shape[0]
    if not all(n_samples == a.shape[0] for a in arrays):
        raise ValueError('the first dimension of all arrays should be equal')
    if indices is None:
        if remove < 1:
            n_remove = round(remove*n_samples)
        else:
            n_remove = remove
        indices = np.random.choice(n_samples, replace=False, size=n_remove)
    return (np.delete(a, indices, axis=0) for a in arrays)

class RobustEstimator(BaseEstimator, ClassifierMixin):
    '''A wrapper to select robust features based on feature reccurence

    Parameters:
    ----------
    estimator: sklearn.base.BaseEstimator object
        Base estimator for feature selection
    
    n_select: int
        Number of features to select

    resample_method: str, choices: ('jackknife', 'bootstrap')
        Resampling method

    remove: float or int
        Fraction (float, 0.0-1.0) or number (int, >= 1) of samples to remove for each round of resampling
    
    max_runs: int
        Maximum rounds of jackknife resampling
    
    recurring_fraction: float
        Minimum required fraction of reccuring samples to select features
    
    grid_search: dict or None
        Parameter grid for optimizing hyper-parameters using GridSearchCV
    
    rfe: bool
        Whether to use RFE to select features rather than select features with higest importance during each run

    Attributes:
    ----------
    feature_recurrence_: array-like, shape (n_features,)
        Number of resampling rounds in which each feature is retained
    
    features_: array-like, shape (n_features,)
        Indices of selected features
    
    feature_selection_matrix_: array-like, (max_runs, n_select), type int
        Indicator matrix for selected features during each run

    feature_importance_matrix_: array-like, (max_runs, n_features)
        Feature importance during each resampling run
    
    feature_importances_mean_: array-like, (n_features,)
        Feature importances averaged across resampling runs
    
    feature_importance_std_: array-like, (n_features,)
        Standard deviation of feature importances across resampling runs
    
    estimator_: sklearn.base.BaseEstimator object
        Estimator fitted on all data using selected features
    
    feature_importances_: array-like, (n_select,)
        Feature importances trained on all data using selected features

    '''
    def __init__(self, estimator, n_select=None, cv=10):
        self.estimator = estimator
        self.n_select = n_select
        self.cv = cv

    def fit(self, X, y, sample_weight=None):
        n_samples, n_features = X.shape
        feature_rank_matrix = np.zeros(n_samples)
        if self.remove == 1:
            max_runs = n_samples
        else:
            max_runs = self.max_runs
        n_select = self.n_select
        if n_select is None:
            n_select = n_features
        feature_rank_matrix = np.zeros((max_runs, n_features))
        feature_importances_matrix = np.zeros((max_runs, n_features))
        # compute sample weight
        if sample_weight is None:
            sample_weight = np.ones(n_samples)
        feature_rank_matrix = np.zeros((max_runs, n_features), dtype=np.int32)
        feature_importances_matrix = np.zeros((max_runs, n_features))

        for train_index, _ in cv.split(X, y, sample_weight):
            estimator = clone(self.estimator)
            estimator.fit(X[train_index], y[train_index], sample_weight=sample_weight[train_index])
            feature_importances = get_feature_importances(estimator)
            feature_importances_matrix.append(feature_importances)
            feature_rank_matrix.append(get_feature_ranking(feature_importances))
            
        estimator = clone(self.estimator)
        estimator.fit(X, y, sample_weight=sample_weight)

        for i_run in range(max_runs):
            if self.resample_method == 'bootstrap':
                X_, y_, sample_weight_ = bootstrap(X, y, sample_weight)
            elif self.resample_method == 'jackknife':  
                indices_remove = None
                if self.remove == 1:
                    indices_remove = i_run
                X_, y_, sample_weight_ = jackknife(X, y, sample_weight,
                    remove=self.remove, indices=indices_remove)
            if self.grid_search is not None:
                cv = GridSearchCV(self.estimator, param_grid=self.grid_search, cv=3)
                cv.fit(X_, y_, sample_weight=sample_weight_)
                self.estimator = cv.best_estimator_
            else:
                self.estimator.fit(X_, y_, sample_weight=sample_weight_)
            if self.rfe:
                rfe = RFE(self.estimator, n_select, step=0.1)
                rfe.fit(X_, y_)
                feature_importances_matrix[i_run] = rfe.support_.astype('float')
                feature_rank_matrix[i_run] = rfe.ranking_
            else:
                if hasattr(self.estimator, 'coef_'):
                    feature_importances = np.square(self.estimator.coef_.flatten())
                else:
                    feature_importances = self.estimator.feature_importances_
                feature_importances_matrix[i_run] = feature_importances
                feature_orders = np.argsort(-feature_importances)
                feature_rank_matrix[i_run, feature_orders] = np.arange(n_features)
        if self.rfe:
            feature_selection_matrix = (feature_rank_matrix == 1).astype(np.int32)
        else:
            feature_selection_matrix = (feature_rank_matrix < n_select).astype(np.int32)
        feature_recurrence = np.sum(feature_selection_matrix, axis=0)
        #robust_features = np.nonzero(feature_recurrence > round(n_samples)*self.recurring_fraction)[0]
        #robust_features = robust_features[np.argsort(feature_recurrence[robust_features])][::-1]
        robust_features = np.argsort(-feature_recurrence)[:n_select]
        feature_selection_matrix = feature_selection_matrix[:, robust_features]

        # refit the estimator
        if self.grid_search is not None:
            cv = GridSearchCV(self.estimator, param_grid=self.grid_search, cv=3)
            cv.fit(X[:, robust_features], y, sample_weight=sample_weight)
            self.estimator = cv.best_estimator_
        else:
            self.estimator.fit(X[:, robust_features], y, sample_weight=sample_weight)

        # save attributes
        self.feature_recurrence_ = feature_recurrence
        self.features_ = robust_features
        self.feature_selection_matrix_ = feature_selection_matrix
        self.feature_importances_matrix_ = feature_importances_matrix
        self.feature_importances_mean_ = np.mean(feature_importances_matrix, axis=0)
        self.feature_importances_std_ = np.std(feature_importances_matrix, axis=0)
        self.estimator_ = self.estimator
        if hasattr(self.estimator, 'coef_'):
            self.feature_importances_ = np.square(self.estimator.coef_.flatten())
        else:
            self.feature_importances_ = self.estimator.feature_importances_
        
        return self

    def predict(self, X):
        return self.estimator_.predict(X)
    
    def score(self, X):
        return self.estimator_.score(X)

class RobustSelector(BaseEstimator, SelectorMixin):
    def __init__(self, estimator, cv=None, n_features_to_select=10, verbose=0):
        self.estimator = estimator
        self.cv = cv
        self.n_features_to_select = n_features_to_select
        self.verbose = verbose
    
    def fit(self, X, y, sample_weight=None):
        n_samples, n_features = X.shape
        # compute sample weight
        if sample_weight is None:
            sample_weight = np.ones(n_samples)
        feature_rank_matrix = []
        feature_importances_matrix = []
        cv = get_splitter(**self.cv)
        for train_index, _ in cv.split(X, y, sample_weight):
            estimator = clone(self.estimator)
            estimator.fit(X[train_index], y[train_index], sample_weight=sample_weight[train_index])
            feature_importances = get_feature_importances(estimator)
            feature_importances_matrix.append(feature_importances)
            feature_rank_matrix.append(get_feature_ranking(feature_importances))
        self.feature_rank_matrix_ = np.vstack(feature_rank_matrix)
        self.feature_importances_matrix_ = np.vstack(feature_importances_matrix)
        self.feature_selection_matrix_ = (self.feature_rank_matrix_ < self.n_features_to_select).astype(np.int32)
        self.feature_recurrence_ = np.mean(self.feature_selection_matrix_, axis=0)

        self.feature_importances_ = self.feature_importances_matrix_.mean(axis=0)
        self.ranking_ = get_feature_ranking(-self.feature_rank_matrix_.mean(axis=0))
        self.support_ = np.zeros(n_features, dtype='bool')
        self.support_[np.argsort(-self.feature_recurrence_)[:self.n_features_to_select]] = True
        return self

    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_

class RpkmFilter(BaseEstimator, SelectorMixin):
    def __init__(self, threshold=1, below=False, pseudo_count=0.001):
        self.threshold = threshold
        self.below = below
        self.pseudo_count = pseudo_count
    
    def set_gene_lengths(self, gene_lengths):
        self.gene_lengths = gene_lengths

    def fit(self, X, y=None):
        if getattr(self, 'gene_lengths') is None:
            raise ValueError('gene_lengths is required for RpkmFilter')
        rpkm = 1e3*X/self.gene_lengths.reshape((1, -1))
        self.rpkm_mean_ = np.exp(np.mean(np.log(rpkm + self.pseudo_count), axis=0)) - self.pseudo_count
        return self
    
    def _get_support_mask(self):
        check_is_fitted(self, 'rpkm_mean_')
        if self.below:
            return self.rpkm_mean_ < self.threshold
        else:
            return self.rpkm_mean_ > self.threshold

class FoldChangeFilter(BaseEstimator, SelectorMixin):
    '''Feature selection based on fold change

    Parameters:
    ----------

    threshold: float
        Fold change threshold
    
    direction: str
        'both': fold change in both direction (up-regulated or down-regulated)
        'up': up-regulated
        'down': down-regulated
    
    below: bool
        True if select features with fold change below threshold
    
    pseudo_count: float
        Pseudo-count to add to input expression matrix

    '''
    def __init__(self, threshold=1, direction='any', below=False, pseudo_count=0.001):
        if threshold <= 0:
            raise ValueError('fold change threshold should be > 0')
        self.direction = direction
        self.threshold = threshold
        self.pseudo_count = pseudo_count
        self.below = below
    
    def fit(self, X, y):
        unique_classes = np.sort(np.unique(y))
        if len(unique_classes) != 2:
            raise ValueError('FoldChangeSelector requires exactly 2 classes, but found {} classes'.format(len(unique_classes)))
        # calculate geometric mean
        X = X + self.pseudo_count
        X_log = np.log(X)
        log_mean = np.zeros(2)
        for i, c in enumerate(unique_classes):
            log_mean[i] = np.mean(X_log[y == c], axis=0)
        logfc = log_mean[1] - log_mean[0]
        if self.direction == 'any':
            self.logfc_ = np.abs(logfc)
        elif self.direction == 'down':
            self.logfc_ = -logfc
        elif self.direction == 'up':
            self.logfc_ = logfc
        else:
            raise ValueError('unknown fold change direction: {}'.format(self.direction))
        return self
    
    def _get_support_mask(self):
        check_is_fitted(self, 'logfc_')
        
        if self.below:
            return self.logfc_ < np.log(self.threshold)
        else:
            return self.logfc_ > np.log(self.threshold)


class ZeroFractionFilter(BaseEstimator, SelectorMixin):
    def __init__(self, max_zero_fraction=0.8):
        self.max_zero_fraction = max_zero_fraction
    
    def fit(self, X, y=None):
        self.zero_fractions_ = np.mean(np.close(X, 0), axis=0)
        return self

    def _get_support_mask(self):
        check_is_fitted(self, 'zero_fractions_')
        return self.zero_fractions_ < self.max_zero_fraction


class LogTransform(BaseEstimator, TransformerMixin):
    def __init__(self, base=None, pseudo_count=0.001):
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
    def __call__(self, estimator, X, y, y_pred, train_index, test_index):
        pass

class CollectMetrics(CVCallback):
    def __init__(self, scoring='roc_auc', classifier='classifier'):
        self.scoring = scoring
        self.metrics = []
        self.classifier = classifier
    
    def __call__(self, estimator, X, y, y_pred, train_index, test_index):
        scorer = get_scorer(self.scoring)
        metrics = {}
        metrics['train_{}'.format(self.scoring)] = scorer(y[train_index], y_pred[train_index])
        if test_index.shape[0] > 1:
            metrics['test_{}'.format(self.scoring)] = scorer(y[test_index], y_pred[test_index])
        self.metrics.append(metrics)
    
    def get_metrics(self):
        if isinstance(self.metrics, list):
            self.metrics = pd.DataFrame.from_records(self.metrics)
        return self.metrics

class CollectPredictions(CVCallback):
    def __init__(self):
        self.predictions = []
    
    def __call__(self, estimator, X, y, y_pred, train_index, test_index):
        self.predictions.append(np.ravel(y_pred))
    
    def get_predictions(self):
        if isinstance(self.predictions, list):
            self.predictions = np.vstack(self.predictions)
        return self.predictions

class FeatureSelectionMatrix(CVCallback):
    def __init__(self, selector='selector'):
        self.matrix = []
        self.selector = selector
    
    def __call__(self, estimator, X, y, y_pred, train_index, test_index):
        if hasattr(estimator, 'selector_'):
            self.matrix.append(estimator.selector_.support_)
    
    def get_matrix(self):
        if isinstance(self.matrix, list):
            self.matrix = np.vstack(self.matrix)
        return self.matrix

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
    classifier_params=None,
    grid_search=None,
    grid_search_param_grid=None,
    grid_search_scoring=None,
    grid_search_cv_params=None
    ):


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
    print('scaler_params: {}'.format(scaler_params))
    if scaler is not None:
        steps.append(('scaler', get_scaler(scaler, **parse_params(scaler_params))))
    if grid_search is not None:
        steps.append(('grid_search', GridSearchCV(
            classifier, param_grid=parse_params(grid_search_), scoring=grid_search_scoring, cv=get_splitter(*grid_search_cv_params))))

    if selector is not None:
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

class CombinedEstimator(BaseEstimator, MetaEstimatorMixin):
    def __init__(self, feature_names=None,
        zero_fraction_filter=None,
        fold_change_filter=None,
        rpkm_filter=None,
        log_transform=None,
        scaler=None,
        scaler_params=None,
        selector=None,
        selector_params=None,
        classifier=None,
        classifier_params=None,
        grid_search=False,
        grid_search_param_grid=None,
        grid_search_scoring=None,
        grid_search_cv_params=None):

        preprocess_steps = []
        if zero_fraction_filter is not None:
            preprocess_steps.append(('zero_fraction_filter', get_selector('zero_fraction_filter', 
                **parse_params(zero_fraction_filter))))
        if rpkm_filter != 'none':
            if feature_names is None:
                raise ValueError('feature_names is required for rpkm_filter')
            gene_lengths = self.get_gene_lengths_from_feature_names(feature_names)
            step = get_selector('rpkm_filter', **parse_params(rpkm_filter))
            step.set_gene_lengths(gene_lengths)
            preprocess_steps.append(('rpkm_filter', step))
        if fold_change_filter != 'none':
            preprocess_steps.append(('fold_change_filter', get_selector(
                'fold_change_filter', **parse_params(fold_change_filter))))
        if log_transform != 'none':
            preprocess_steps.append(('log_transform', get_scaler(
                'log_transform', **parse_params(log_transform))))
        if scaler is not None:
            preprocess_steps.append(('scaler', get_scaler(scaler, **parse_params(scaler_params))))
        self.preprocess_steps = preprocess_steps

        self.grid_search = grid_search
        self.grid_search_param_grid = parse_params(grid_search_param_grid)
        self.grid_search_scoring = grid_search_scoring
        self.grid_search_cv_params = parse_params(grid_search_cv_params)
        self.classifier_name = classifier
        self.classifier_params = parse_params(classifier_params)
        self.selector_name = selector
        self.selector_params = parse_params(selector_params)
    
    @staticmethod
    def get_gene_lengths_from_feature_names(feature_names):
        feature_info = pd.Series(feature_names).str.split('|', expand=True)
        feature_info.columns = ['gene_id', 'gene_type', 'gene_name', 'feature_id', 'transcript_id', 'start', 'end']
        feature_info['start'] = feature_info['start'].astype('int')
        feature_info['end'] = feature_info['end'].astype('int')
        feature_info['length'] = feature_info['end'] - feature_info['start']
        return feature_info['length'].values
    
    def fit(self, X, y=None, sample_weight=None):
        X_new = X
        self.features_ = np.arange(X.shape[1])
        for name, step in self.preprocess_steps:
            X_new = step.fit_transform(X_new, y)
            setattr(self, name + '_', step)
            if isinstance(step, SelectorMixin):
                self.features_ = self.features_[step.get_support()]
        self.classifier_ = get_classifier(self.classifier_name, **self.classifier_params)
        if self.grid_search:
            self.grid_search_ = GridSearchCV(estimator=self.classifier_, 
                param_grid=self.grid_search_param_grid, 
                scoring=self.grid_search_scoring, 
                cv=get_splitter(**self.grid_search_cv_params))
            self.grid_search_.fit(X_new, y)
            self.classfier_ = self.grid_search_.best_estimator_
            self.best_classifier_params_ = self.grid_search_.best_params_
            self.classifier_.set_params(**self.grid_search_.best_params_)
        if self.selector_name:
            self.selector_ = get_selector(self.selector_name, estimator=self.classifier_, **self.selector_params)
            X_new = self.selector_.fit_transform(X_new, y)
            self.features_ = self.features_[self.selector_.get_support()]
        self.classifier_.fit(X_new, y)
        return self
    
    def transform(self, X, y=None):
        check_is_fitted(self, 'classifier_')
        for name, step in self.preprocess_steps:
            X = step.transform(X)
        if self.selector_name:
            X = self.selector_.transform(X)
        return X
    
    def predict(self, X):
        X = self.transform(X)
        return self.classifier_.predict(X)

def cross_validation(estimator, X, y, params=None, sample_weight=None, callbacks=None):
    splitter = get_splitter(**params)
    logger.info('params: {}'.format(params))
    logger.info(str(splitter))
    logger.info(str(splitter.get_n_splits(X, y)))
    pbar = tqdm(unit='split', total=splitter.get_n_splits())
    for index in splitter.split(X, y):
        train_index, test_index = index
        estimator.fit(X[train_index], y[train_index])
        y_pred = estimator.predict(X)
        for callback in callbacks:
            callback(estimator, X, y, y_pred, train_index, test_index)
        pbar.update(1)
    pbar.close()