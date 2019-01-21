# Frequently Asked Questions

## How to install and set up the environment
We have packaged everything into docker and singularity. The easiest and recommended way is via Singularity.

You can refer to this doc's [Singularity part](https://github.com/lulab/exSeek/blob/master/docs/requirements.md)

## How to use exSeek
exSeek is an integrative tool for exRNA processing and feature selection. We use snakemake for parallel running and further integrate snakemake pipeline into one single command.

Details of preparing steps are described [here](https://github.com/lulab/exSeek/blob/master/README.md). Basically you should complete the following steps before running the command:

* Install exseek and requirements
* Prepare genome and annotation
* prepare input files in right file path
* set up configuration

Then you can run the command, you can specify the module you want to run and dataset you provide.

```
${exseek_path}/bin/exseek.py quality_control --dataset ${dataset}
```

## What is Snakemake
The [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system is a tool to create reproducible and scalable data analyses. We have hide the details of snakemake and you only need to run one single command. However you can customize some of the codes if you are familiar with snakemake.

## How to set configurations in config file
There are many parameters to be specified. You should make a new copy of config file in config directory. For example you can nake one copy of scirep.yaml. Then rename the file to config/${dataset}.yaml. 

Other parameters are defined in snakemake/default_config.yaml. You may also change parameters. 

## How to generate report
After running some modules, e.g., mapping, normalization and evaluation. You can open jupyter notebook files in notebooks file. The only thing to do is to fill in the dataset name and sequencing type.

For example:

```
dataset='scirep'
sequencing_type = 'short'
```

Then you can get plots of your mapping, processing and feature selection details in one jupyter notebook.

**Note:** the notebook is based on exseek output style. If you process your data on your own without exseek and only need the jupyter to generate plots, you should change the codes for file paths in jupyter notebook to successfully generate plots.


## When bugs appear:
The quickest way is to create a [new issue](https://github.com/lulab/exSeek/issues) 

If you want us to add more functions in exseek, please create a [new issue](https://github.com/lulab/exSeek/issues)

## Why did some jobs occasionally fail with no extra error message except 'CalledProcessError'?

The most possible cause is no available memory. You can confirm the problem by running the linux command `dmesg` or open
and examine the most recent message. Message like "Out of memory: Kill process ** or sacrifice child" clearly indicates that memory problem occurred.
Some jobs (e.g. mapping using STAR, samtools sort, bedtools sort)requires large amount memory especially
 when the input number of reads is large. Try to reduce the number of parallels
with the `-j` option or set memory limit in `config/cluster.yaml` for a particular job if you run the jobs on a
computer cluster. 

## How to rerun downstream steps from a specific step in a pipeline

Sometimes we need to rerun a pipeline from a step, usually after changing the configuration file.
Snakemake is not aware of changes in configuration file and we need to rerun the pipeline by ourselves.
The `--forcerun` option in snakemake allows rerunning a step and all steps that depend on the output files
of the step. For example, to rerun the `count_matrix` step, just run:

```bash
exseek.py count_matrix -d $dataset --forcerun count_matrix
```

