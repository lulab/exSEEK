# File format

## Input files

**`config/${dataset}.yaml`**

A full list of parameters and descriptions can be found in the default configuration file [default_config.yaml](../config/default_config.yaml).
This file is writtern in [YAML](https://yaml.org) format.

**`data/${dataset}/sample_classes.txt`**

Tab-separated text file with two columns: sample_id, sample_class. The first line should be columns names: sample_id, sample_class.

**`data/compare_groups.yaml`**

This YAML file contains definition of sample groups for feature selection. In the following example, there are two comparisons: "Normal-Cancer", "Normal-stage_A".
In "Normal-Cancer", the negative class and positive class are "Normal" and samples with "stage_A", "stage_B", "stage_C" respectively. 
In "Normal-stage_A", "Normal" samples are compared to only "stage_A" samples. the class labels after comparison name should match class labels defined in `data/${dataset}/sample_classes.txt`. 

```yaml
Normal-Cancer: ['Normal', 'stage_A,stage_B,stage_C']
Normal-stage_A: ['Normal', 'stage_A']
```

**`data/batch_info.txt`**

Tab-separated text file with at least two columns: sample_id, batch1, batch2. The first line should be columns names.
Note that only one batch variable is supported at one time, although it is possible to specify multiple batch variables in the file.
When multiple batch variables are specified, you can set the `batch_index` variable in `config/${dataset}.yaml` file.

**`output/${dataset}/count_matrix/count_matrix.txt`**

Tab-separated text file of read counts. Columns are samples and rows are features. The first line should be column names. The first column should be row names.


## Feature selection

**Variables in file patterns**

| Variable | Descrpition |
| :--- | :--- |
| `output_dir` | Output directory for the dataset, e.g. `output/dataset` |
| `preprocess_method` | Combination of matrix processing methods |
| `count_method` | Type of feature counts, e.g. `domains_combined`, `domains_long`, `transcript`, `featurecounts` |
| `compare_group` | Name of the negative-positive class pair defined in `compare_groups.yaml` |
| `classifier` | Classifier defined in the configuration file |
| `n_select` | Maximum number of features to select |
| `selector` | Feature selection method, e.g. `robust`, `rfe` |
| `fold_change_filter_direction` | Direction of fold change for filtering features. Three possible values: `up`, `down` and `any` |

**List of files in output directory**

| File name pattern | Descrpition |
| :--- | :--- |
| `features.txt` | Selected features. Plain text with one column: feature names |
| `feature_importances.txt` | Plain text with two columns: feature name, feature importance |
| `samples.txt` | Sample IDs in input matrix selected for feature selection |
| `classes.txt` | Sample class labels selected for feature selection |
| `final_model.pkl` | Final model fitted on all samples in Python pickle format |
| `metrics.train.txt` | Evaluation metrics on training data. First row is metric names. First column is index of each train-test split |
| `metrics.test.txt` | Same format with `metrics.train.txt` on test data. |
| `cross_validation.h5` | Cross-validation details in HDF5 format. |

**Cross validation details \(cross\_validation.h5\)**

| Dataset name | Dimension | Description |
| :--- | :--- | :--- |
| feature\_selection | \(n\_splits, n\_features\) | Binary matrix indicating features selected in each cross-validation split |
| labels | \(n\_samples,\) | True class labels |
| predicted\_labels | \(n\_splits, n\_samples\) | Predicted class labels on all samples |
| predictions | \(n\_splits, n\_samples\) | Predicted probabilities of the positive class \(or decision function for SVM\) |
| train\_index | \(n\_splits, n\_samples\) | Binary matrix indicating training samples in each cross-validation split |
