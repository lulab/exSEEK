# exSeek
[![Build Status](https://travis-ci.com/lulab/exSeek-dev.svg?token=CyRgUWsqWCctKvAxMXto&branch=master)](https://travis-ci.com/lulab/exSeek-dev)

## Workflow

![workflow](assets/wholepipe.png)

## Installation

Install required software packages according to [requirements](docs/requirements.md)

## Prepare genome and annotations

### Download and process genome sequences

Refer to the [documentation](docs/genome_and_annotations.md) for details.

### Extract GTFs and generate mapping indexes
```bash
snakemake --snakefile snakemake/prepare_genome.snakemake \
    --configfile snakemake/config.yaml \
    --rerun-incomplete -k
```

## Input files

## Configuration

All parameters are specified in a configuration file in [YAML](https://en.wikipedia.org/wiki/YAML) format.

An example configuration file is `snakemake/config.yaml`.

The parameter values in the configuration file can also be overrided through the `--config` option in [snakemake](https://snakemake.readthedocs.io/en/stable/executable.html).

The following parameters should be changed:

| Parameter | Description | Example |
| ------ | ----------- | ------- |
| `genome_dir` | Directory for genome and annotation files | `genome/hg38` |
| `data_dir` | Directory for input files | `data/scirep` |
| `temp_dir` | Temporary directory | `tmp` |
| `sample_id_file` | A text file containing sample IDs | `data/scirep/sample_ids.txt` |
| `output_dir` | Directory for all output files | `output/scirep` |
| `tools_dir` | Directory for third-party tools |
| `aligner` | Mapping software | `bowtie2` |
| `adaptor` | 3' adaptor sequence for single-end RNA-seq | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` |
| `python2` | Path to Python 2 | `/apps/anaconda2/bin/python` |
| `python3` | Path to Python 3 | `/apps/anaconda2/bin/python` |

### Command line arguments for snakemake

| Option | Description |
| ------ | ----------- |
| --config | Additional configuration parameters |
| -j | Number of parallel jobs |
| --dryrun | Do not execute |
| -k | Do not stop when an independent job fails |

## Mapping (small RNA-seq)

### Generate snakemake rules for sequential mapping
```bash
bin/generate_snakemake.py sequential_mapping --rna-types rRNA,miRNA,piRNA,Y_RNA,srpRNA,tRNA,snRNA,snoRNA,lncRNA,mRNA,tucpRNA \
    -o snakemake/mapping_small/sequential_mapping.snakemake
```

```bash
snakemake --snakefile snakemake/mapping_small.snakemake \
    --configfile snakemake/config.yaml \
    --rerun-incomplete -k
```

**Output files**

| File name | Descrpition |
| --------- | ----------- |
| `output/${dataset}/cutadapt/${sample_id}.fastq` | Reads with adaptor trimmed |
| `output/${dataset}/tbam/${sample_id}/${rna_type}.bam` | BAM files in transcript coordinates |
| `output/${dataset}/gbam/${sample_id}/${rna_type}.bam` | BAM files in genome coordinates |
| `output/${dataset}/unmapped/${sample_id}/${rna_type}.fa.gz` | Unmapped reads in each step |
| `output/${dataset}/fastqc/${sample_id}_fastqc.html` | FastQC report file |
| `output/${dataset}/summary/fastqc.html` | Summary report for FastQC |
| `output/${dataset}/summary/read_counts.txt` | Summary table for read counts |


## Generate expression matrix
```bash
snakemake --snakefile snakemake/expression_matrix.snakemake \
    --configfile snakemake/config.yaml \
    --rerun-incomplete -k
```

## Call domains for long RNA
```bash
snakemake --snakefile snakemake/call_domains_long.snakemake \
    --configfile snakemake/config.yaml \
    --rerun-incomplete -k
```