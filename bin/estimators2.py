import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin, MetaEstimatorMixin, is_classifier
from sklearn.base import clone
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.svm import LinearSVC
from sklearn.metrics import roc_auc_score, accuracy_score, roc_curve, \
    precision_score, recall_score, f1_score, \
    average_precision_score
from sklearn.feature_selection.base import SelectorMixin
from sklearn.preprocessing import StandardScaler, RobustScaler, MinMaxScaler, MaxAbsScaler
from sklearn.model_selection import KFold, StratifiedKFold, ShuffleSplit, LeaveOneOut, \
        RepeatedKFold, RepeatedStratifiedKFold, LeaveOneOut, StratifiedShuffleSplit
from sklearn.utils.class_weight import compute_sample_weight
from sklearn.model_selection import GridSearchCV, check_cv
from sklearn.feature_selection import RFE, RFECV, SelectFromModel
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_X_y
from abc import ABC, ABCMeta, abstractmethod
from tqdm import tqdm
import pickle
import json
import os
import subprocess
import shutil
from scipy.cluster import hierarchy
from scipy.spatial.distance import pdist, squareform
from tqdm import tqdm
import logging
import copy
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
        except json.JSONDecodeError as e:
            pass
        finally:
            return params
    else:
        return {}

def search_dict(data, keys):
    '''Intersect given keys with a dict and return the subset

    Parameters:
        data: dict
        keys: list of keys to search
    Returns:
        dict
    '''
    return {key:data[key] for key in keys if key in data}

def predict_proba(estimator, X):
    try:
        proba = estimator.predict_proba(X)
    except AttributeError:
        proba = estimator.decision_function(X)
    return proba

def get_scorer(scoring):
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
    if name == 'logistic_regression':
        return LogisticRegression(**search_dict(params, 
            ('penalty', 'dual', 'C', 'tol', 'fit_intercept', 'solver',
                'class_weight', 'max_iter', 'n_jobs', 'random_state', 'verbose')))
    elif name == 'logistic_regression_l1':
        return LogisticRegression(penalty='l1', **search_dict(params, 
            ('dual', 'C', 'tol', 'fit_intercept', 'solver',
            'class_weight', 'max_iter', 'n_jobs', 'random_state', 'verbose')))
    elif name == 'logistic_regression_l2':
        return LogisticRegression(penalty='l2', **search_dict(params, 
            ('dual', 'C', 'tol', 'fit_intercept', 'solver',
            'class_weight', 'max_iter', 'n_jobs', 'random_state', 'verbose')))
    elif name == 'random_forest':
        return RandomForestClassifier(**search_dict(params,
            ('n_estimators', 'criterion', 'max_depth', 'min_samples_split', 'min_samples_leaf', 
             'min_weight_fraction_leaf', 'max_features', 'max_leaf_nodes', 
             'min_impurity_decrease', 'min_impurity_split', 'oob_score',
             'n_jobs', 'verbose', 'random_state', 'class_weight')))
    elif name == 'linear_svm':
        return LinearSVC(**search_dict(params,
            ('penalty', 'loss', 'dual', 'tol', 'C', 'fit_intercept', 
             'intercept_scaling', 'class_weight', 'verbose',
             'random_state', 'max_iter')))
    elif name == 'decision_tree':
        return DecisionTreeClassifier(**search_dict(params,
            ('criterion', 'splitter', 'max_depth', 'min_samples_split', 'min_samples_leaf',
             'min_weight_fraction_leaf', 'max_features', 'max_leaf_nodes', 'min_impurity_decrease',
             'min_impurity_split')))
    elif name == 'extra_trees':
        return ExtraTreesClassifier(**search_dict(params,
            ('n_estimators', 'criterion', 'max_depth', 'min_samples_split', 'min_samples_leaf', 
             'min_weight_fraction_leaf', 'max_features', 'max_leaf_nodes', 
             'min_impurity_decrease', 'min_impurity_split', 'oob_score',
             'n_jobs', 'verbose', 'random_state', 'class_weight')))
    else:
        raise ValueError('unknown classifier: {}'.format(name))

def get_selector(name, estimator=None, n_features_to_select=None, **params):
    if name == 'robust':
        return RobustSelector(estimator, n_features_to_select=n_features_to_select, **search_dict(params,
         ('cv', 'verbose')))
    elif name == 'max_features':
        return SelectFromModel(estimator, threshold=-np.inf, max_features=n_features_to_select)
    elif name == 'feature_importance_threshold':
        return SelectFromModel(estimator, **search_dict(params, 'threshold'))
    elif name == 'rfe':
        return RFE(estimator, n_features_to_select=n_features_to_select, **search_dict(params, 
        ('step', 'verbose')))
    elif name == 'rfecv':
        return RFECV(estimator, n_features_to_select=n_features_to_select, **search_dict(params,
         ('step', 'cv', 'verbose')))
    elif name == 'fold_change_filter':
        return FoldChangeFilter(**search_dict(params,
        ('threshold', 'direction', 'below', 'pseudo_count')))
    elif name == 'zero_fraction_filter':
        return ZeroFractionFilter(**search_dict(params,
        ('threshold',)))
    elif name == 'rpkm_filter':
        return RpkmFilter(**search_dict(params,
        ('threshold',)))
    elif name == 'rpm_filter':
        return RpmFilter(**search_dict(params,
        ('threshold',)))
    elif name == 'diffexp_filter':
        return DiffExpFilter(**search_dict(params,
        ('threshold', 'max_features', 'script', 'temp_dir', 'score_type', 'diffexp_method')))
    elif name == 'null':
        return NullSelector()
    else:
        raise ValueError('unknown selector: {}'.format(name))

def get_splitter(random_state=None, **params):
    splitter = params.get('splitter')
    if splitter is None:
        return check_cv(**params)
    if splitter == 'kfold':
        return KFold(random_state=random_state, **search_dict(params, ('n_splits', 'shuffle')))
    elif splitter == 'stratified_kfold':
        return StratifiedKFold(random_state=random_state, **search_dict(params, ('n_splits', 'shuffle')))
    elif splitter == 'repeated_stratified_kfold':
        return RepeatedStratifiedKFold(random_state=random_state, **search_dict(params, ('n_splits', 'n_repeats')))
    elif splitter == 'shuffle_split':
        return ShuffleSplit(random_state=random_state, **search_dict(params, ('n_splits', 'test_size', 'train_size')))
    elif splitter == 'stratified_shuffle_split':
        return StratifiedShuffleSplit(random_state=random_state, **search_dict(params, ('n_splits', 'test_size', 'train_size')))
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

class RobustSelector(BaseEstimator, SelectorMixin):
    '''Feature selection based on recurrence

    Parameters:
    ----------

    estimator: object
        A classifier that provides feature importances through feature_importances_ or coef_ attribute after calling the fit method.
    
    cv: int or splitter object
        Specifies how to subsample the original dataset
    
    n_features_to_select: int
        Maximum number of features to select
    
    Attributes:
    ----------
    
    support_: array-like, shape (n_features,)
        Boolean mask indicating features selected
    
    ranking_: array-like, shape (n_features,)
        Ranking of feature importances starting from 0 to n_features - 1.
        Smaller ranks indicates higher importance.
    
    feature_recurrence_: array-like, shape (n_features,)
        Number of times each feature is selected across resampling runs divided by total number of resampling runs.
    
    feature_selection_matrix_: array-like, shape (n_splits, n_features)
        A boolean matrix indicates features selected in each resampling run
    
    feature_rank_matrix_: array-like, shape (n_splits, n_features)
        Feature ranks in each resampling run
    
    feature_importances_matrix_: array-like, shape (n_splits, n_features)
        Feature importances in each resampling run
    
    feature_importances_: array-like, shape (n_features,)
        Average feature importances across resampling runs
    '''
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
    def __init__(self, threshold=1, below=False, pseudo_count=0.01):
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
        if self.below:
            self.support_ = self.rpkm_mean_ < self.threshold
        else:
            self.support_ = self.rpkm_mean_ > self.threshold
        return self
    
    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_
        
class RpmFilter(BaseEstimator, SelectorMixin):
    '''Feature selection based on geometric mean expression value across samples (RPM)

    Parameters:
    ----------

    threshold: float
        Expression value threshold
    
    below: bool
        True if select features with expression value below threshold
    
    pseudo_count: float
        Pseudo-count to add to input expression matrix during calculating the geometric mean
    
    Attributes:
    ----------
    
    support_: bool | array-like, shape (n_features,)
        Boolean mask indicating features selected

    '''
    def __init__(self, threshold=1, below=False, pseudo_count=0.01):
        self.threshold = threshold
        self.below = below
        self.pseudo_count = pseudo_count

    def fit(self, X, y=None):
        self.rpm_mean_ = np.exp(np.mean(np.log(X + self.pseudo_count), axis=0)) - self.pseudo_count
        if self.below:
            self.support_ = self.rpm_mean_ < self.threshold
        else:
            self.support_ = self.rpm_mean_ > self.threshold
        return self
    
    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_

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
    def __init__(self, threshold=1, direction='any', below=False, pseudo_count=0.01):
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
        X_log = np.log2(X)
        log_mean = np.zeros((2, X.shape[1]))
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
        
        if self.below:
            self.support_ = self.logfc_ < np.log(self.threshold)
        else:
            self.support_ = self.logfc_ > np.log(self.threshold)
        return self
    
    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_

class ZeroFractionFilter(BaseEstimator, SelectorMixin):
    '''Feature selection based on fraction of zero values

    Parameters:
    ----------

    threshold: float
        Features with zero values above this fraction will be filtered out
    
    eps: float
        Define zero values as values below this number

    '''
    def __init__(self, threshold=0.8, eps=0.0):
        self.threshold = threshold
        self.eps = eps
    
    def fit(self, X, y=None):
        self.zero_fractions_ = np.mean(X <= self.eps, axis=0)
        self.support_ = self.zero_fractions_ < self.threshold
        return self

    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_

class HvgFilter(BaseEstimator, SelectorMixin):
    def __init__(self, threshold=0):
        self.threshold = threshold
    
    def fit(self, X, y=None):
        mean = np.mean(X, axis=0)
        std = np.std(X, axis=0)

class FisherDiscriminantRatioFilter(BaseEstimator, SelectorMixin):
    '''Feature selection based on fraction of zero values

    Parameters:
    ----------

    threshold: float
        Features with zero values above this fraction will be filtered out
    
    eps: float
        Define zero values as values below this number
    
    Attributes:
    ----------

    support_: array-like, shape (n_features,)
        Boolean mask indicating features selected
    
    scores_: array-like, shape (n_features,)
        Fisher's discriminant ratios
    '''
    def __init__(self, threshold=0):
        self.threshold = threshold
    
    def fit(self, X, y):
        unique_classes = np.unique(y)
        if len(unique_classes) != 2:
            raise ValueError('Fisher discriminant ratio requires 2 classes')
        unique_classes = np.sort(unique_classses)
        mean = np.zeros(2)
        var = np.zeros(2)
        for i in range(2):
            mean[i] = np.mean(X[y == unique_classes[i]], axis=0, ddof=1)
            var[i] = np.var(X[y == unique_classes[i]], axis=0, ddof=1)
        self.scores_ = (mean[1] - mean[0])/(var[0] + var[1])
        self.support_ = self.scores_ > self.threshold
    
    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_

class DiffExpFilter(BaseEstimator, SelectorMixin):
    '''Feature selection based on differential expression

    Parameters:
    ----------

    threshold: float
        Features with zero values above this fraction will be filtered out
    
    max_features: int
        Maximum number of features to select
    
    score_type: str
        Scores for ranking features.
        Allowed values: 'neglogp', 'logfc', 'pi_score'
        'pi_score'[2]: | log FC - log (padj) |
    
    method: str
        Differential expression method to use.
        Available methods: 'deseq2', 'edger_glmlrt', 'edger_glmqlf', 'edger_exact', 'wilcox'.
    
    temp_dir: str
        Temporary directory for storing input and output files for differential expression
    
    script: str
        Path of the script (to differential_expression.R)
        The script takes two input files: matrix.txt and sample_classes.txt and outputs a table named results.txt.
        The output file contains at least two columns: log2FoldChange, padj
    
    fold_change_direction: str
        Direction of fold change filter
        Allowed values: up, down or any.
    
    fold_change_threshold: float
        Threshold for absolute fold change
    
    References:
    ----------

    1. Rosario, S.R., Long, M.D., Affronti, H.C., Rowsam, A.M., Eng, K.H., and Smiraglia, D.J. (2018). 
        Pan-cancer analysis of transcriptional metabolic dysregulation using The Cancer Genome Atlas. Nature Communications 9, 5330.
    2. Xiao, Y., Hsiao, T.-H., Suresh, U., Chen, H.-I.H., Wu, X., Wolf, S.E., and Chen, Y. (2014). 
        A novel significance score for gene selection and ranking. Bioinformatics 30, 801–807.

    
    Attributes:
    ----------

    support_: array-like, shape (n_features,)
        Boolean mask indicating features selected
    
    padj_: array-like, shape (n_features,)
        Adjusted p-values for each feature
    
    logfc_: array-like, shape (n_features,)
        Log2 fold change

    '''
    def __init__(self, threshold=0, max_features=None, 
            score_type='adjusted_pvalue', temp_dir=None,
            script=None, method='deseq2',
            fold_change_direction='any',
            fold_change_threshold=1):
        self.threshold = threshold
        self.max_features = max_features
        self.score_type = score_type
        self.temp_dir = temp_dir
        self.method = method
        self.script = script
        self.fold_change_direction = fold_change_direction
        self.fold_change_threshold = fold_change_threshold
    
    def fit(self, X, y):
        if self.temp_dir is None:
            raise ValueError('parameter temp_dir is required for DiffExpSelector')
        if self.script is None:
            raise ValueError('parameter script is required for DiffExpSelector')
        # save expression matrix to file
        matrix = pd.DataFrame(X.T, 
            index=['F%d'%i for i in range(X.shape[1])],
            columns=['S%d'%i for i in range(X.shape[0])])
        if not os.path.isdir(self.temp_dir):
            logger.debug('create diffexp dir: {}'.format(self.temp_dir))
            os.makedirs(self.temp_dir)
        try:
            matrix_file = os.path.join(self.temp_dir, 'matrix.txt')
            logger.debug('write expression matrix to : {}'.format(matrix_file))
            matrix.to_csv(matrix_file, sep='\t', na_rep='NA', index=True, header=True)
            # save sample_classes to file
            sample_classes = pd.Series(y, index=matrix.columns.values)
            sample_classes = sample_classes.map({0: 'negative', 1: 'positive'})
            sample_classes.name = 'label'
            sample_classes.index.name = 'sample_id'
            sample_classes_file = os.path.join(self.temp_dir, 'sample_classes.txt')
            logger.debug('write sample classes to: {}'.format(sample_classes_file))
            sample_classes.to_csv(sample_classes_file, sep='\t', na_rep='NA', index=True, header=True)
            output_file = os.path.join(self.temp_dir, 'results.txt')
            logger.debug('run differential expression script: {}'.format(self.script))
            subprocess.check_call([self.script, '--matrix', matrix_file, 
                '--classes', sample_classes_file,
                '--method', self.method,
                '--positive-class', 'positive', '--negative-class', 'negative',
                '-o', output_file])
            # read results
            logger.debug('read differential expression results: {}'.format(output_file))
            results = pd.read_table(output_file, sep='\t', index_col=0)
        finally:
            # remove temp_dir after reading the output file or an exception occurs
            logger.debug('remove diffexp directory: {}'.format(self.temp_dir))
            shutil.rmtree(self.temp_dir, ignore_errors=True)
        self.logfc_ = results.loc[:, 'log2FoldChange']
        self.padj_ = results.loc[:, 'padj']
        if self.score_type == 'neglogp':
            self.scores_ = -np.log10(self.padj_)
        elif self.score_type == 'logfc':
            self.scores_ = np.abs(self.logfc_)
        elif self.score_type == 'pi_value':
            self.scores_ = np.abs(self.logfc_)*(-np.log10(self.padj_))
        else:
            raise ValueError('unknown differential expression score type: {}'.format(self.score_type))
        # apply fold change filter
        fc_support = np.ones(X.shape[1], dtype='bool')
        if self.fold_change_direction != 'any':
            logger.debug('apply fold change filter: {} > {:.2f}'.format(self.fold_change_direction, self.fold_change_threshold))
            if self.fold_change_threshold == 'up':
                fc_support = self.logfc_ > np.log2(self.fold_change_threshold)
            elif self.fold_change_threshold == 'down':
                fc_support = self.logfc_ < -np.log2(self.fold_change_threshold)
            else:
                raise ValueError('unknown fold change direction: {}'.format(self.fold_change_direction))
        # compute support mask
        if self.max_features is not None:
            # sort feature scores in descending order and get top features
            indices = np.nonzero(fc_support)[0]
            indices = indices[np.argsort(-self.scores_[indices])][:self.max_features]
            self.support_ = np.zeros(X.shape[1], dtype='bool')
            self.support_[indices] = True
        else:
            # select features with scores above a given threshold
            self.support_ = (self.scores_ > self.threshold) & fc_support
        return self
    
    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_

class NullSelector(BaseEstimator, SelectorMixin):
    '''A null selector that select all features

    Attributes:
    ----------

    support_: array-like, shape (n_features,)
        Boolean mask indicating features selected
    '''
    def fit(self, X, y=None):
        self.support_ = np.ones(X.shape[1], dtype='bool')
        return self
    
    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_
    

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
        scorer = get_scorer(self.scoring)
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
        if hasattr(estimator, 'selector_'):
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

class CombinedEstimator(BaseEstimator, MetaEstimatorMixin):
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
        grid_search_params=None):
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

        X_new = X
        self.features_ = np.arange(X.shape[1])
        for name, step in self.preprocess_steps_:
            X_new = step.fit_transform(X_new, y)
            setattr(self, name + '_', step)
            if isinstance(step, SelectorMixin):
                self.features_ = self.features_[step.get_support()]
        logger.debug('add classifier "{}" with parameters: {}'.format(self.classifier, self.classifier_params))
        self.classifier_ = get_classifier(self.classifier, **self.classifier_params)
        if self.grid_search:
            logger.debug('add grid_search with parameters: {}'.format(self.grid_search_params))
            grid_search_params = copy.deepcopy(self.grid_search_params)
            if 'cv' in grid_search_params:
                grid_search_params['cv'] = get_splitter(**grid_search_params['cv'])
            grid_search_params['param_grid'] = grid_search_params['param_grid'][self.classifier]
            self.grid_search_ = GridSearchCV(estimator=self.classifier_, 
                **search_dict(grid_search_params, ('param_grid', 'scoring', 'cv', 
                'fit_params', 'verbose', 'return_train_score', 'error_score')))
            self.grid_search_.fit(X_new, y)
            self.classfier_ = self.grid_search_.best_estimator_
            self.best_classifier_params_ = self.grid_search_.best_params_
            self.classifier_.set_params(**self.grid_search_.best_params_)
        if self.selector:
            logger.debug('add selector "{}" with parameters: {}'.format(self.selector, self.selector_params))
            logger.debug('number of features to select: {}'.format(self.n_features_to_select))
            self.selector_ = get_selector(self.selector, estimator=self.classifier_, 
                n_features_to_select=self.n_features_to_select, **self.selector_params)
            X_new = self.selector_.fit_transform(X_new, y)
            self.features_ = self.features_[self.selector_.get_support()]
        self.classifier_.fit(X_new, y)
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

def cross_validation(estimator, X, y, sample_weight='auto', params=None, callbacks=None):
    splitter = get_splitter(**params)
    logger.debug('start cross-validation')
    logger.debug('cross-validation parameters: {}'.format(params))
    logger.debug('number of cross-validation splits: {}'.format(splitter.get_n_splits(X, y)))
    pbar = tqdm(unit='split', total=splitter.get_n_splits(X, y))
    for index in splitter.split(X, y):
        train_index, test_index = index
        estimator = clone(estimator)
        sample_weight_ = sample_weight
        if sample_weight == 'auto':
            sample_weight_ = compute_sample_weight(class_weight=None, y=y[train_index])
        estimator.fit(X[train_index], y[train_index], sample_weight=sample_weight_)
        y_pred_labels = estimator.predict(X)
        y_pred_probs = estimator.predict_proba(X)
        for callback in callbacks:
            callback(estimator, X, y, y_pred_labels, y_pred_probs, train_index, test_index)
        pbar.update(1)
    pbar.close()
