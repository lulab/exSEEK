# Feature Selection

## Run feature selection module

Assume that we have already run `normalization` module and selected best matrix processing method based on the UCA score, 
we can run feature selection module using the following command:
```bash
exseek.py feature_selection -d ${dataset}
```

## Output files

### Outuput directory

Feature selection results using one combination of parameters are saved in a separate directory:
```
${output_dir}/cross_validation/${preprocess_method}.${count_method}/${compare_group}/${classifier}.${n_select}.${selector}.${fold_change_filter_direction}
```

**Variables in file patterns**

| Variable | Descrpition |
| --------- | ----------- |
| `output_dir` | Output directory for the dataset, e.g. `output/dataset` |
| `preprocess_method` | Combination of matrix processing methods |
| `count_method` | Type of feature counts, e.g. `domains_combined`, `domains_long`, `transcript`, `featurecounts` |
| `compare_group` | Name of the negative-positive class pair defined in `compare_groups.yaml` |
| `classifier` | Classifier defined in the configuration file |
| `n_select` | Maximum number of features to select |
| `selector` | Feature selection method, e.g. `robust`, `rfe` |
| `fold_change_filter_direction` | Direction of fold change for filtering features. Three possible values: `up`, `down` and `any` |

### Files in output directory

| File name pattern | Descrpition |
| --------- | ----------- |
| `features.txt` | Selected features. Plain text with one column: feature names |
| `feature_importances.txt` | Plain text with two columns: feature name, feature importance |
| `samples.txt` | Sample IDs in input matrix selected for feature selection |
| `classes.txt` | Sample class labels selected for feature selection |
| `final_model.pkl` | Final model fitted on all samples in Python pickle format |
| `metrics.train.txt` | Evaluation metrics on training data. First row is metric names. First column is index of each train-test split |
| `metrics.test.txt` | Same format with `metrics.train.txt` on test data. |
| `cross_validation.h5` | Cross-validation details in HDF5 format. |

**Cross validation details (cross_validation.h5)**

| Dataset name | Dimension | Description |
| ------------ | --------- | ----------- |
| feature_selection | (n_splits, n_features) | Binary matrix indicating features selected in each cross-validation split |
| labels | (n_samples,) | True class labels |
| predicted_labels | (n_splits, n_samples) | Predicted class labels on all samples |
| predictions | (n_splits, n_samples) | Predicted probabilities of the positive class (or decision function for SVM) |
| train_index | (n_splits, n_samples) | Binary matrix indicating training samples in each cross-validation split |
