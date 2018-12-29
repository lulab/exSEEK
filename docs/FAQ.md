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

### What is Snakemake
The [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow management system is a tool to create reproducible and scalable data analyses. We have hide the details of snakemake and you only need to run one single command. However you can customize some of the codes if you are familiar with snakemake.

### How to set configurations in config file
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


