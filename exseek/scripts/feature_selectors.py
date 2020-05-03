from sklearn.base import BaseEstimator, ClassifierMixin, TransformerMixin, MetaEstimatorMixin, is_classifier
from sklearn.base import clone
from sklearn.feature_selection.base import SelectorMixin
from sklearn.model_selection import GridSearchCV, check_cv
from sklearn.feature_selection import RFE, RFECV, SelectFromModel
from sklearn.neural_network import MLPClassifier
import numpy as np
import pandas as pd
from sklearn.utils.validation import check_is_fitted
from sklearn.utils import check_X_y
from tqdm import tqdm
import os
import subprocess
import shutil
from copy import deepcopy
import logging
from utils import search_dict, get_feature_importances, python_args_to_r_args
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger(__name__)

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

class MaxFeaturesSelector(BaseEstimator, SelectorMixin):
    '''Select given number of features from model

    Parameters:
    ----------

    estimator: object
        A classifier that provides feature importances through feature_importances_ or coef_ attribute after calling the fit method.
    
    grid_search: bool
        Whether to optimize hyper-parameters of the estimator by grid search
    
    grid_search_params: dict
        Parameters passed to GridSearchCV
    
    n_features_to_select: int
        Maximum number of features to select
    
    Attributes:
    ----------

    support_: array-like, shape (n_features,)
        Boolean mask indicating features selected
    
    feature_importances_: array-like, shape (n_features,)
        Average feature importances across resampling runs
    '''
    def __init__(self, estimator, n_features_to_select=10, 
        grid_search=False, grid_search_params=None):
        self.estimator = estimator
        self.n_features_to_select = n_features_to_select
        self.grid_search = grid_search
        self.grid_search_params = grid_search_params
    
    def fit(self, X, y, sample_weight=None):
        if self.grid_search is not None:
            grid_search = GridSearchCV(self.estimator,
                **self.grid_search_params)
            grid_search.fit(X, y, sample_weight=sample_weight)
            self.estimator_ = grid_search.best_estimator_
            self.best_classifier_params_ = grid_search.best_params_
            self.estimator_.set_params(**self.best_classifier_params_)
        self.estimator_.fit(X, y, sample_weight=sample_weight)
        self.feature_importances_ = get_feature_importances(self.estimator)
        self.support_ = np.zeros(X.shape[1], dtype='bool')
        self.support_[np.argsort(-self.feature_importances_)][:self.n_features_to_select] = True
    
    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_
    
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
        for train_index, _ in self.cv.split(X, y, sample_weight):
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
        self.support_[np.argsort(-self.feature_recurrence_)][:self.n_features_to_select] = True
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

    def fit(self, X, y=None, **kwargs):
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

    def fit(self, X, y=None, **kwargs):
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
    
    def fit(self, X, y, **kwargs):
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
    
    def fit(self, X, y=None, **kwargs):
        self.zero_fractions_ = np.mean(X <= self.eps, axis=0)
        self.support_ = self.zero_fractions_ < self.threshold
        return self

    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_

class HvgFilter(BaseEstimator, SelectorMixin):
    def __init__(self, threshold=0):
        self.threshold = threshold
    
    def fit(self, X, y=None, **kwargs):
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
    
    def fit(self, X, y, **kwargs):
        unique_classes = np.unique(y)
        if len(unique_classes) != 2:
            raise ValueError('Fisher discriminant ratio requires 2 classes')
        unique_classes = np.sort(unique_classes)
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
    def __init__(self, 
            rscript_path='Rscript',
            threshold=0, 
            max_features=None, 
            score_type='adjusted_pvalue',
            temp_dir=None,
            script=None, method='deseq2',
            fold_change_direction='any',
            fold_change_threshold=1):
        self.rscript_path = rscript_path
        self.threshold = threshold
        self.max_features = max_features
        self.score_type = score_type
        self.temp_dir = temp_dir
        self.method = method
        self.script = script
        self.fold_change_direction = fold_change_direction
        self.fold_change_threshold = fold_change_threshold
    
    def fit(self, X, y, **kwargs):
        if self.temp_dir is None:
            raise ValueError('parameter temp_dir is required for DiffExpFilter')
        if self.script is None:
            raise ValueError('parameter script is required for DiffExpFilter')
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
            subprocess.check_call([self.rscript_path, self.script, 
                '--matrix', matrix_file, 
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

class SIS(BaseEstimator, SelectorMixin):
    '''Sure Independence Screening

    Original R function:

    SIS(x, y, family = c("gaussian", "binomial", "poisson", "cox"),
        penalty = c("SCAD", "MCP", "lasso"), concavity.parameter = switch(penalty,
        SCAD = 3.7, 3), tune = c("bic", "ebic", "aic", "cv"), nfolds = 10,
        type.measure = c("deviance", "class", "auc", "mse", "mae"),
        gamma.ebic = 1, nsis = NULL, iter = TRUE, iter.max = ifelse(greedy ==
        FALSE, 10, floor(nrow(x)/log(nrow(x)))), varISIS = c("vanilla", "aggr",
        "cons"), perm = FALSE, q = 1, greedy = FALSE, greedy.size = 1,
        seed = 0, standardize = TRUE)

    Parameters:
    ----------

    temp_dir: str
        directory for storing temporary files
    
    sis_args: dict
        Arguments to pass to SIS function
    

    '''
    def __init__(self, rscript_path='Rscript', temp_dir=None, n_features_to_select=None, sis_params=None):
        self.rscript_path = rscript_path
        self.temp_dir = temp_dir
        self.n_features_to_select = n_features_to_select
        self.sis_params = sis_params
        if n_features_to_select is not None:
            self.sis_params['nsis'] = n_features_to_select
        if self.sis_params is None:
            self.sis_params.update(deepcopy(sis_params))
        if self.temp_dir is None:
            raise ValueError('temp_dir is required for SIS')
    
    def fit(self, X, y):
        if not os.path.isdir(self.temp_dir):
            os.makedirs(self.temp_dir)
        try:
            pd.DataFrame(X).to_csv(os.path.join(self.temp_dir, 'matrix.txt'), 
                sep='\t', header=False, index=False, na_rep='NA')
            pd.Series(y).to_csv(os.path.join(self.temp_dir, 'labels.txt'), header=False, index=False)
            #print(y[:10])
            r_script = r'''
library(SIS)
X <- read.table('{temp_dir}/matrix.txt', header=FALSE, check.names=FALSE)
X <- as.matrix(X)
y <- read.table('{temp_dir}/labels.txt')[,1]
y <- as.numeric(y)
model <- SIS(X, y, {sis_params})
write.table(model$coef.est, '{temp_dir}/coef.txt', sep='\t', col.names=FALSE, row.names=TRUE, quote=FALSE)
write.table(model$ix, '{temp_dir}/ix.txt', col.names=FALSE, row.names=FALSE, quote=FALSE)
    '''
            print(self.sis_params)
            r_script = r_script.format(temp_dir=self.temp_dir, sis_params=python_args_to_r_args(self.sis_params))
            print(python_args_to_r_args(self.sis_params))
            script_file = os.path.join(self.temp_dir, 'run_SIS.R')
            with open(script_file, 'w') as f:
                f.write(r_script)
            logger.debug('execute R script: ' + r_script)
            subprocess.check_call([self.rscript_path, script_file], shell=False)
            # read outputs
            #coef = pd.read_table(os.path.join(self.temp_dir, 'coef.txt'), sep='\t', index=True, header=None)
            # read indices of selected features
            ix = pd.read_table(os.path.join(self.temp_dir, 'ix.txt'), sep='\t', header=None)
            indices = ix.iloc[:, 0].values - 1
            self.support_ = np.zeros(X.shape[1], dtype='bool')
            self.support_[indices] = True
        finally:
            pass
            #shutil.rmtree(self.temp_dir)
    
    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_

class RandomSubsetSelector(BaseEstimator, SelectorMixin):
    '''Large scale feature selection based on max feature weights on feature subsets

    Parameters:
    ----------

    estimator: BaseEstimator object
        Internal estimator to use for feature selection
    
    n_subsets: int
        Number of random feature subsets
    
    subset_size: int
        Number of features in each subset
    
    n_features_to_select: int
        Maximum number of features to select
    
    random_state: RandomState
        Random number generator
    '''
    def __init__(self, estimator, n_subsets=40, subset_size=50, n_features_to_select=10, random_state=None):
        self.estimator = estimator
        self.n_subsets = n_subsets
        self.subset_size = subset_size
        self.n_features_to_select = n_features_to_select
        self.random_state = random_state
        
    def fit(self, X, y, sample_weight=None):
        n_features = X.shape[1]
        feature_weights = np.zeros((self.n_subsets, n_features))
        rng = np.random.RandomState(self.random_state)
        for subset_index in range(self.n_subsets):
            subset = rng.choice(n_features, size=self.subset_size, replace=False)
            estimator = clone(self.estimator)
            estimator.fit(X[:, subset], y, sample_weight=sample_weight)
            feature_weights[subset_index, subset] = get_feature_importances(estimator)
        # get local maximum
        feature_weights = np.max(feature_weights, axis=0)
        self.features_ = np.argsort(-feature_weights)[:self.n_features_to_select]
        self.support_ = np.zeros(n_features, dtype='bool')
        self.support_[self.features_] = True
    
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
    def fit(self, X, y=None, **kwargs):
        self.support_ = np.ones(X.shape[1], dtype='bool')
        return self
    
    def _get_support_mask(self):
        check_is_fitted(self, 'support_')
        return self.support_
    
def get_selector(name, estimator=None, n_features_to_select=None, **params):
    if name == 'RobustSelector':
        return RobustSelector(estimator, n_features_to_select=n_features_to_select, **search_dict(params,
         ('cv', 'verbose')))
    elif name == 'MaxFeatures':
        return SelectFromModel(estimator, threshold=-np.inf, max_features=n_features_to_select)
    elif name == 'RandomSubsetSelector':
        return RandomSubsetSelector(estimator, n_features_to_select=n_features_to_select, **search_dict(params,
        ('n_subsets', 'subset_size', 'random_state')))
    elif name == 'FeatureImportanceThreshold':
        return SelectFromModel(estimator, **search_dict(params, 'threshold'))
    elif name == 'RFE':
        return RFE(estimator, n_features_to_select=n_features_to_select, **search_dict(params, 
        ('step', 'verbose')))
    elif name == 'RFECV':
        return RFECV(estimator, n_features_to_select=n_features_to_select, **search_dict(params,
         ('step', 'cv', 'verbose')))
    elif name == 'FoldChangeFilter':
        return FoldChangeFilter(**search_dict(params,
        ('threshold', 'direction', 'below', 'pseudo_count')))
    elif name == 'ZeroFractionFilter':
        return ZeroFractionFilter(**search_dict(params,
        ('threshold',)))
    elif name == 'RpkmFilter':
        return RpkmFilter(**search_dict(params,
        ('threshold',)))
    elif name == 'RpmFilter':
        return RpmFilter(**search_dict(params,
        ('threshold',)))
    elif name == 'DiffExpFilter':
        return DiffExpFilter(max_features=n_features_to_select, **search_dict(params,
        ('rscript_path', 'threshold', 'script', 'temp_dir', 'score_type', 'method')))
    elif name == 'ReliefF':
        from skrebate import ReliefF
        return ReliefF(n_features_to_select=n_features_to_select,
            **search_dict(params, ('n_jobs', 'n_neighbors', 'discrete_limit')))
    elif name == 'SURF':
        from skrebate import SURF
        return SURF(n_features_to_select=n_features_to_select,
            **search_dict(params, ('n_jobs', 'discrete_limit')))
    elif name == 'MultiSURF':
        from skrebate import MultiSURF
        return MultiSURF(n_features_to_select=n_features_to_select,
            **search_dict(params, ('n_jobs', 'discrete_limit')))
    elif name == 'SIS':
        return SIS(n_features_to_select=n_features_to_select, 
            **search_dict(params, ('rscript_path', 'temp_dir', 'sis_params')))
    elif name == 'NullSelector':
        return NullSelector()
    else:
        raise ValueError('unknown selector: {}'.format(name))