# exSeek
[![Build Status](https://travis-ci.com/lulab/exSeek-dev.svg?token=CyRgUWsqWCctKvAxMXto&branch=master)](https://travis-ci.com/lulab/exSeek-dev)

## Workflow

![workflow](assets/whole_pipe.png)

## Installation

Install required software packages according to [requirements](docs/requirements.md)

Download the scripts:

```bash
git clone https://github.com/lulab/exSeek-dev.git
```

## Prepare genome and annotations

**Processed files**: `/BioII/lulab_b/shared/genomes/hg38`

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
| `${input_dir}/fastq/${sample_id}_1.fastq`, `${input_dir}/fastq/${sample_id}_2.fastq` | Read files (paired-end sequencing) |
| `${input_dir}/sample_ids.txt` | A text file with one sample ID per line. |
| `${input_dir}/sample_classes.txt` | A tab-deliminated file (with header) with two columns: sample_id, label |
| `${input_dir}/batch_info.txt` | A comma-deliminated file (with header) with at least two columns: sample_id, batch1, batch2, ... |
| `${input_dir}/reference_genes.txt` | A text file with reference gene IDs. |
| `${input_dir}/compare_groups.yaml` | A YAML file defining positive and negative classes. |

**compare_groups.yaml**

Every key-value pairs defines a compare group and a negative-positive class pair:
```yaml
Normal-CRC: ["Healthy Control", "Colorectal Cancer"]
```

## Configuration

All parameters are specified in a configuration file in [YAML](https://en.wikipedia.org/wiki/YAML) format.

An example configuration file is (snakemake/config.yaml).

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

## Submit jobs to a computer cluster using snakemake

Please refer the [link](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html#cluster-configuration) for descriptions of cluster configuration file.

### IBM LSF

**Configuration file**: `snakemake/cluster.yaml`

Here is an example configuration:

```yaml
__default__:
  queue: Z-LU
  name: {rule}.{wildcards}
  stderr: logs/cluster/{rule}/{wildcards}.stderr
  stdout: logs/cluster/{rule}/{wildcards}.stdout
  threads: {threads}
  resources: span[hosts=1]
```

**Commonly used parameters**

| Parameter | Description |
| ------ | ----------- |
| `__default__` | Rule name (`__default__`) for default configuration) | 
| `queue` | Queue name |
| `name` | Job name |
| `stderr` | Log file for standard error |
| `stdout` | Log file for standard output |
| `threads` | Number of parallel threads for a job |
| `resources` | Resource requirements. `span[hosts=1]` prevents parallel jobs from being submitted to different nodes |

**Run snakemake**
```bash
snakemake --snakefile snakemake/${snakefile} \
    --configfile snakemake/config.yaml \
    --cluster 'bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} \
-o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}' \
    --cluster-config snakemake/cluster.yaml \
    --rerun-incomplete -k -j40
```
**Note**: replace `${snakefile}` with a Snakefile.

## Quality control

```bash
snakemake --snakefile snakemake/quality_control.snakemake \
    --configfile snakemake/config.yaml \
    --rerun-incomplete -k
```

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
| `${output_dir}/cutadapt/${sample_id}.fastq` | Reads with adaptor trimmed |
| `${output_dir}/tbam/${sample_id}/${rna_type}.bam` | BAM files in transcript coordinates |
| `${output_dir}/gbam/${sample_id}/${rna_type}.bam` | BAM files in genome coordinates |
| `${output_dir}/unmapped/${sample_id}/${rna_type}.fa.gz` | Unmapped reads in each step |
| `${output_dir}/fastqc/${sample_id}_fastqc.html` | FastQC report file |
| `${output_dir}/summary/fastqc.html` | Summary report for FastQC (HTML) |
| `${output_dir}/summary/fastqc.txt`  | Summary table for FastQC |
| `${output_dir}/summary/fastqc.ipynb` | Summary report for FastQC (Jupyter notebook) |
| `${output_dir}/summary/read_counts.txt` | Summary table for read counts |
| `${output_dir}/stats/mapped_read_length_by_sample/${sample_id}` | Length distribution of mapped reads |


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

**Count matrix**

* File path: `${output_dir}/count_matrix/transcript.txt`
* First row: sample IDs
* First column: feature names
* Feature name: `gene_id|gene_type|gene_name`


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


**Read count matrix**

* File path: `${output_dir}/count_matrix/domain_long.txt`
* First row: sample IDs
* First column: feature names
* Feature name: `gene_id|gene_type|gene_name|domain_id|transcript_id|start|end`

## Normalization

### Output files

| File name | Description |
| `${output_dir}/normalized_matrix/${normalization_method}.${imputation_method}.${batch_removal_method}.txt` |
| `${output_dir}/matrix_processing/normalization/${normalization_method}.txt` |
| `${output_dir}/matrix_processing/imputation/${normalization_method}.${imputation_method}.txt` |
| `${output_dir}/matrix_processing/batch_removal/${batch_removal_method}.${batch_index}.txt` |
