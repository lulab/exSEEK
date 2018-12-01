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

| File name | Description |
| ------ | ----------- |
| `${input_dir}/fastq/${sample_id}.fastq` | Read files (single-end sequencing) |
| `${input_dir}/fastq/${sample_id}_1.fastq`, `${input_dir}/fastq/${sample_id}_1.fastq` | Read files (paired-end sequencing) |
| `${input_dir}/sample_ids.txt` | A text file with one sample ID per line. |

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

### Output files

| File name | Descrpition |
| --------- | ----------- |
| `snakemake/sequential_mapping.snakemake` | Snakefile for sequential mapping. Required by snakemake/mapping_small.snakemake |
| `output/${dataset}/cutadapt/${sample_id}.fastq` | Reads with adaptor trimmed |
| `output/${dataset}/tbam/${sample_id}/${rna_type}.bam` | BAM files in transcript coordinates |
| `output/${dataset}/gbam/${sample_id}/${rna_type}.bam` | BAM files in genome coordinates |
| `output/${dataset}/unmapped/${sample_id}/${rna_type}.fa.gz` | Unmapped reads in each step |
| `output/${dataset}/fastqc/${sample_id}_fastqc.html` | FastQC report file |
| `output/${dataset}/summary/fastqc.html` | Summary report for FastQC (HTML) |
| `output/${dataset}/summary/fastqc.txt`  | Summary table for FastQC |
| `output/${dataset}/summary/fastqc.ipynb` | Summary report for FastQC (Jupyter notebook) |
| `output/${dataset}/summary/read_counts.txt` | Summary table for read counts |
| `output/${dataset}/stats/mapped_read_length_by_sample/${sample_id}` | Length distribution of mapped reads |


## Generate expression matrix
```bash
snakemake --snakefile snakemake/expression_matrix.snakemake \
    --configfile snakemake/config.yaml \
    --rerun-incomplete -k
```

### Output files
| File name | Descrpition |
| --------- | ----------- |
| `${output_dir}/count_matrix/transcript.txt` | Count matrix of transcripts |
| `${output_dir}/count_matrix/htseq.txt` | Count matrix of genes generated using HTSeq-count |
| `${output_dir}/count_matrix/featurecounts.txt` | Count matrix of genes generated using featureCounts |
| `${output_dir}/counts_by_biotype/${count_method}/${sample_id}/${rna_type}` | Gene/transcript counts generated using a feature counting tool |



## Call domains for long RNA

```bash
snakemake --snakefile snakemake/call_domains_long.snakemake \
    --configfile snakemake/config.yaml \
    --rerun-incomplete -k
```

### Output files

| File name | Descrpition |
| --------- | ----------- |
| `${output_dir}/domain_counts/${bin_size}/${pvalue}/${sample_id}.bed` | Read counts in long RNA domains (BED format with read counts in Column 5 |
| `${output_dir}/count_matrix/domain_${pvalue}.txt` | Read count matrix of long RNA domains |
| `${output_dir}/domains/${bin_size}/${pvalue}.bed` | Long RNA domain locations |
| `${output_dir}/domains_recurrence/${bin_size}/${pvalue}.bed` | Recurrence of long RNA domains among samples (Column 5) |

