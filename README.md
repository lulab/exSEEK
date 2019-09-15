# exSEEK

**exRNA Biomarker Discovery for Liquid Biopsy**

> **Note:**
> * The exSEEK program starts from a data matrix of gene expression (read counts of each gene in each sample) and performs normalization, feature selection and evaluation. 
> * Meanwhile, we provide some pipelines and QC steps for the [pre-process](https://github.com/lulab/exSEEK_docs/tree/master/pre-process) of exRNA-seq (including long and short  cfRNA-seq/exoRNA-seq) raw data. 
> * We also recommend other alternatives for the pre-process, such as [exceRpt](https://github.com/gersteinlab/exceRpt), that is specifically developed for the process of exRNA-seq raw reads.


**Table of Contents:**

* [Installation](#istallation)
* [Usage](#Usage)
  * [Input files](#input-files)
  * [Normalization](#normalization)
  * [Feature selection](#feature-selection)
  * [Advanced usage](#advanced-usage)
* [Copyright and License Information](#copyright-and-license-information)
* [Citation](#citation)

---


## Installation

For easy installation, you can use the [exSEEK image](https://hub.docker.com/r/ltbyshi/exseek) of [docker](https://www.docker.com) with all dependencies installed:

```bash
docker pull ltbyshi/exseek
```

Alternatively, you can use use [singularity](https://singularity.lbl.gov/) or [udocker](https://github.com/indigo-dc/udocker) to run the container for Linux kernel < 3 or if you don't have permission to use docker.

## Usage

Run the main program `exseek.py` from docker:

```bash
docker run --rm -it -v $PWD:/workspace -w /workspace ltbyshi/exseek exseek.py
```

The exSEEK directory was cloned to /apps/exseek in the docker.

You can create a bash script named `exseek` and set the script executable:
```bash
#! /bin/bash
docker run --rm -it -v $PWD:/workspace -w /workspace ltbyshi/exseek exseek.py "$@"
```

After adding the file to one of the directory in the `$PATH` variable, you can simply run: `exseek`.


The basic usage of exSEEK is:

```bash
exseek ${step_name} -d ${dataset}
```

> **Note:**
> * Other arguments are passed to [snakemake](https://snakemake.readthedocs.io/en/stable/)
> * Specify number of processes to run in parallel with *-j*
> * `${step_name}` is one of `normalization` and `cross_validation`.
> * `${dataset}` is the name of your dataset that should match the prefix of your configuration file described in the following section.


### Input files

An example can be found in `example_data` directory with the following structure:
```
example_data/
├── config
│   └── example.yaml
├── data
│   └── example
│       ├── batch_info.txt
│       ├── compare_groups.yaml
│       ├── sample_classes.txt
│       └── sample_ids.txt
└── output
    └── example
        └── count_matrix
            └── mirna_and_domains_rna.txt
```

> **Note:**
> * `config/example.yaml`: configuration file
> * `data/example/batch_info.txt`: table of batch information
> * `data/example/compare_groups.yaml`: configuration file for definition of positive and negative samples
> * `data/example/sample_classes.txt`: table of sample labels
> * `output/example/count_matrix/mirna_and_domains_rna.txt`: input matrix of read counts

You can create your own data directory with the above directory structure. 
Multiple datasets can be put in the same directory by replacing "example" with your own dataset names.

More information about input and output files can be found on [File Format](docs/file_format.md) page.

### Normalization

Run:

```bash
exseek normalization -d ${dataset}
```

This will generate normalized expression matrix for every combination of methods with the following file name pattern:

`output/${dataset}/matrix_processing/filter.${imputation_method}.Norm_${normalization_method}.Batch_${batch_removal_method}_${batch_index}.${count_method}.txt`

You can specify normalization methods by setting the value of `normalization_method` and the batch removal methods
by setting the value of `batch_removal_method` in in `config/${dataset}.yaml`.

Supported normalization methods: TMM, RLE, CPM, CPM_top, UQ, null

Supported batch removal methods: limma, ComBat, RUV, null

When the method name is set to "null", the step is skipped.

`${batch_index}` is the column number (start from 1) in `config/${dataset}/batch_info.txt` to be used to remove batch effects.

### Feature selection

Run:

```bash
exseek feature_selection -d ${dataset}
```

This will evaluate all combinations of feature selection methods and classifiers by cross-validation.

Three summary files will be generated:

* `output/${dataset}/summary/cross_validation/metrics.test.txt`
* `output/${dataset}/summary/cross_validation/metrics.train.txt`
* `output/${dataset}/summary/cross_validation/feature_stability.txt`

Cross-validation results and trained models for individual combinations are in this directory:

`output/${dataset}/feature_selection/filter.${imputation_method}.Norm_${normalization_method}.Batch_${batch_removal_method}_${batch_index}.${count_method}/${compare_group}/${classifier}.${n_select}.${selector}.${fold_change_filter_direction}`

Selected list of features are in `features.txt`.

> **Note:**
> More information about output files can be found on [File format](docs/file_format.md) page. Detailed parameters of feature selection and classifiers can be found in [config/machine_learning.yaml](config/machine_learning.yaml).


### Advanced Usage

* [Click here](docs/README.md) to see details

## Copyright and License Information

Copyright (C) 2019 Tsinghua University, Beijing, China 

This program is licensed with commercial restriction use license. Please see the [LICENSE](https://github.com/lulab/exSEEK_docs/blob/master/LICENSE) file for details.

## Citation

Binbin Shi, Jingyi Cao, Xupeng Chen and Zhi John Lu (2019) exSEEK: an integrative computational framework for identifying extracellular RNA biomarkers in liquid biopsy

