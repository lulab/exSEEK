# Installation

## Install through docker (recommended)

exSEEK requires a lot of other tools to run. To save efforts for environment setup, you can use our prepared docker image:
To get more information about docker, you can visit the [official site](https://www.docker.com) of docker.

After installing docker, you can download the image through the following command:

```bash
docker pull ltbyshi/exseek
```

Then you can invoke the main script `exseek` from docker:

```bash
docker run --rm -it -v $PWD:/workspace -w /workspace ltbyshi/exseek exseek.py "$@"
```

This will mount the current directory to `/workspace` as root directory of input and output files. You can also replace `$PWD` with another directory 

Alternatively, you can use use [singularity](https://singularity.lbl.gov/) or [udocker](https://github.com/indigo-dc/udocker) to run the container for Linux kernel < 3 or if you don't have permission to use docker.

### Install through singularity

First install singularity according to (https://singularity.lbl.gov/install-linux).

Then build a singularity image from docker:

```bash
mkdir -p singularity
singularity build singularity/exseek.simg docker://ltbyshi/exseek
```

Finally you need to change the variable `container.image` to the path of the singularity image in the configuration file.

To run a pipeline through udocker, you can add `--runner singularity option to the `exseek` command.

### Install through udocker

First install udocker according to (https://github.com/indigo-dc/udocker/blob/master/doc/installation_manual.md).

Then download the docker image using udocker:

```bash
udocker pull ltbyshi/exseek
udocker create --name=exseek ltbyshi/exseek
```

Finally you need to change the variable `container.image` to the path of the name of the udocker image in the configuration file.

To run a pipeline through udocker, you can add `--runner udocker` option to the `exseek` command.

## Manual installation

### Install required software

Before running exSEEK, you need to install the following software. It is recommended to install most of the software through [conda](https://docs.conda.io/en/latest/) 
with the [bioconda](https://bioconda.github.io/) channel and [conda-forge](https://conda-forge.org/) channel appended.

You can refer to the [Dockerfile](../docker/Dockerfile) for commands to install the dependencies. You should also ensure that the executables are added to the `$PATH` environment variable.


| Software  | Version | Description  |
|-----------|---------|--------------|
| Python  | 3.6  | Python interpreter  |
| R  | 3.5.3  | R interpreter  |
| Java |  | Java runtime support |
| pandas  |   | Python package for dataframe operations.  |
| matplotlib |  | Python package for plotting. |
| seaborn |  | Python package for high-level plotting. |
| h5py |  | Python package for matrix storage. |
| scikit-learn |  | Python package for machine learning |
| mlextend |  | Python package for extra machine learning algorithms compatible with scikit-learn |
| skrebate |  | Python package for feature selection compatible with scikit-learn |
| flask |  | Python package for building lightweight web apps |
| jinja2 |  | Python package for template rendering |
| umap |  | Python package for dimensional reduction using UMAP |
| snakemake |  | Python package for pipeline building |
| tqdm |  | Python package for progress monitoring |
| bedtools |  | Tools for operation of BED files. |
| samtools |  | Tools for operation of SAM/BAM files. |
| STAR |  | RNA mapping software. |
| bowtie2 |  | RNA mapping software. |
| subread |  | Read counting (featureCounts) |
| RSEM |  | Gene quantification |
| bamtools |  | Tools for operation of BAM files. |
| cutadapt |  | Versatile tool for adapter trimming. |
| picard |  | Tools for operation of SAM/BAM files. |
| gffread |  | Tool for conversion between GTF/GFF and other formats. |
| gffcompare |  | Tool for merging GTF/GFF files. |
| bedToGenePred |  | Convert BED file to GenePred format |
| genePredToGtf |  | Convert genePred file to GTF format |
| bedGraphToBigWig |  | Convert bedGraph file to BigWig format |
| bigWigToBedGraph |  | Convert BigWig file to BedGraph format |
| HT-seq |  | Read counting |
| FastX Toolkit |  | Tools for operation of FASTQ/FASTA files |
| BioPython |  | Python package for reading/writing various file formats used in bioinformatics |
| MultiQC |  | Python package for generating graphical reports for quality control results |
| pigz |  | Parallel version of gzip |

### Install exSEEK

The next command installs exSEEK as a Python package, with a main script named `exseek` added to PATH:
```bash
python setup.py install
```

After installation, you can test by running:

```bash
exseek --help
```
