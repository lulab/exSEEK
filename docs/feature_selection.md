# Robust Feature Selection 


## Command line usage
```
usage: feature_selection.py evaluate [-h] --matrix MATRIX --sample-classes
                                     SAMPLE_CLASSES [--transpose]
                                     [--positive-class POSITIVE_CLASS]
                                     [--negative-class NEGATIVE_CLASS]
                                     [--use-log]
                                     [--scaler {zscore,robust,max_abs,min_max,none}]
                                     [--method {logistic_regression,random_forest,linear_svm}]
                                     --output-dir OUTPUT_DIR
                                     [--compute-sample-weight]
                                     [--splitter {kfold,stratified_kfold,shuffle_split,repeated_stratified_kfold,leave_one_out,stratified_shuffle_split}]
                                     [--rfe] [--rfe-step RFE_STEP]
                                     [--n-select N_SELECT]
                                     [--n-splits N_SPLITS]
                                     [--n-repeats N_REPEATS]
                                     [--test-size TEST_SIZE] [--scorer SCORER]
                                     [--robust-select]
                                     [--robust-resample-method {bootstrap,jackknife}]
                                     [--robust-max-runs ROBUST_MAX_RUNS]
                                     [--robust-jackknife-remove ROBUST_JACKKNIFE_REMOVE]
                                     [--top-features-by-median TOP_FEATURES_BY_MEDIAN]
                                     [--remove-zero-features REMOVE_ZERO_FEATURES]

optional arguments:
  -h, --help            show this help message and exit
  --matrix MATRIX, -i MATRIX
                        input feature matrix (rows are samples and columns are
                        features
  --sample-classes SAMPLE_CLASSES
                        input file containing sample classes with 2 columns:
                        sample_id, sample_class
  --transpose           transpose the feature matrix
  --positive-class POSITIVE_CLASS
                        comma-separated list of sample classes to use as
                        positive class
  --negative-class NEGATIVE_CLASS
                        comma-separates list of sample classes to use as
                        negative class
  --use-log             apply log2 to feature matrix
  --scaler {zscore,robust,max_abs,min_max,none}
                        method for scaling features
  --method {logistic_regression,random_forest,linear_svm}
                        feature selection method
  --output-dir OUTPUT_DIR, -o OUTPUT_DIR
                        output directory
  --compute-sample-weight
                        compute sample weight to balance classes
  --splitter {kfold,stratified_kfold,shuffle_split,repeated_stratified_kfold,leave_one_out,stratified_shuffle_split}
  --rfe                 use RFE to select features
  --rfe-step RFE_STEP   number/fraction of features to eliminate in each step
  --n-select N_SELECT   number of features to select
  --n-splits N_SPLITS   number of splits for kfold, stratified_kfold and
                        shuffle_splits
  --n-repeats N_REPEATS
                        number of repeats for repeated_stratified_kfold and
                        repeated_kfold
  --test-size TEST_SIZE
                        fraction/number of samples for testing
  --scorer SCORER       metric to use
  --robust-select       use robust feature selection by selecting recurring
                        features across resampling runs
  --robust-resample-method {bootstrap,jackknife}
                        resampling method for robust feature selection
  --robust-max-runs ROBUST_MAX_RUNS
                        number of resampling runs for robust feature
                        selections
  --robust-jackknife-remove ROBUST_JACKKNIFE_REMOVE
                        number/fraction of samples to remove during each
                        resampling run for robust feature selection
  --top-features-by-median TOP_FEATURES_BY_MEDIAN
                        select this number of features with highest median
                        values across samples
  --remove-zero-features REMOVE_ZERO_FEATURES
                        remove features that have fraction of zero values
                        above this value
```


