# exSEEK

> exRNA Biomarker Discovery for Liquid Biopsy
>
> Date: "4.22.2019"

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Contact](#contact)
- [Copyright and License Information](#copyright-and-license-information)


## Installation

Install required software packages according to [requirements](docs/requirements.md)

Download the scripts:

```bash
git clone https://github.com/lulab/exSeek-dev.git
```


## Usage

Run `exseek.py --help` to get basic usage:

```
usage: exseek.py [-h] --dataset DATASET [--config-dir CONFIG_DIR] [--cluster]
                 [--cluster-config CLUSTER_CONFIG]
                 [--cluster-command CLUSTER_COMMAND]
                 [--singularity SINGULARITY]
                 [--singularity-wrapper-dir SINGULARITY_WRAPPER_DIR]
                 {quality_control,prepare_genome,mapping,count_matrix,call_domains,normalization,feature_selection,update_sequential_mapping,update_singularity_wrappers}

exSeek main program

positional arguments:
  {quality_control,prepare_genome,mapping,count_matrix,call_domains,normalization,feature_selection,update_sequential_mapping,update_singularity_wrappers}

optional arguments:
  -h, --help            show this help message and exit
  --dataset DATASET, -d DATASET
                        dataset name
  --config-dir CONFIG_DIR, -c CONFIG_DIR
                        directory for configuration files
  --cluster             submit to cluster
  --cluster-config CLUSTER_CONFIG
                        cluster configuration file ({config_dir}/cluster.yaml
                        by default)
  --cluster-command CLUSTER_COMMAND
                        command for submitting job to cluster (default read
                        from {config_dir}/cluster_command.txt
  --singularity SINGULARITY
                        singularity image file
  --singularity-wrapper-dir SINGULARITY_WRAPPER_DIR
                        directory for singularity wrappers
```

> **Note**
> * Other arguments are passed to *snakemake*
> * Specify number of processes to run in parallel with *-j*


## Contact

[GitHub Issues](https://github.com/lulab/exSEEK/issues)

Binbin Shi <>


## Copyright and License Information
Copyright (C) 2019 Tsinghua University of California, Beijing, China 

Authors: Binbin Shi, Xupeng Chen, ..., and Zhi John Lu 

This program is licensed with commercial restriction use license. Please see the attached LICENSE file for details.
