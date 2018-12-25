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
| `${input_dir}/sample_classes.txt` | A tab-deliminated file (with header) with two columns: sample_id, label (optional) |
| `${input_dir}/batch_info.txt` | A comma-deliminated file (with header) with at least two columns: sample_id, batch1, batch2, ... (optional) |
| `${input_dir}/compare_groups.yaml` | A YAML file defining positive and negative classes. (optional) |
| `${config_dir}/${dataset}.yaml` | A YAML file for configuration parameters for the dataset |

**compare_groups.yaml**

Every key-value pairs defines a compare group and a negative-positive class pair:
```yaml
Normal-CRC: ["Healthy Control", "Colorectal Cancer"]
```

## Configuration

All parameters are specified in a configuration file in [YAML](https://en.wikipedia.org/wiki/YAML) format.

The default configuration file is (snakemake/default_config.yaml).

The parameter values in the configuration file can also be overrided through the `--config` option in [snakemake](https://snakemake.readthedocs.io/en/stable/executable.html).

The following parameters should be changed:

| Parameter | Description | Example |
| ------ | ----------- | ------- |
| `genome_dir` | Directory for genome and annotation files | `genome/hg38` |
| `data_dir` | Directory for input files | `data/scirep` |
| `temp_dir` | Temporary directory | `tmp` |
| `output_dir` | Directory for all output files | `output/scirep` |
| `tools_dir` | Directory for third-party tools |
| `aligner` | Mapping software | `bowtie2` |
| `adaptor` | 3' adaptor sequence for single-end RNA-seq | `AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` |

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

**cluster configuration**: `config/cluster.yaml`

Here is an example configuration:

```yaml
__default__:
  queue: queue
  name: {rule}.{wildcards}
  stderr: logs/cluster/{rule}/{wildcards}.stderr
  stdout: logs/cluster/{rule}/{wildcards}.stdout
  threads: {threads}
  resources: span[hosts=1]
```

**cluster command**: `config/cluster_command.txt`

```
bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}
```

**Commonly used parameters**

| Parameter | Description |
| ------ | ----------- |
| `__default__` | Rule name (`__default__`) for default configuration) | 
| `queue` | Queue name (required) |
| `name` | Job name |
| `stderr` | Log file for standard error |
| `stdout` | Log file for standard output |
| `threads` | Number of parallel threads for a job |
| `resources` | Resource requirements. `span[hosts=1]` prevents parallel jobs from being submitted to different nodes |

Refer to the snakemake [documentation](https://snakemake.readthedocs.io/en/stable/).

**Run snakemake**
```bash
${exseek_path}/bin/exseek.py --step ${step} --dataset ${dataset} --cluster -j40
```
**Note**: replace `${snakefile}` with a Snakefile.

## Quality control

```bash
${exseek_path}/bin/exseek.py --step quality_control --dataset ${dataset}
```

## Mapping (small RNA-seq)

### Generate snakemake rules for sequential mapping
```bash
bin/generate_snakemake.py sequential_mapping --rna-types rRNA,miRNA,piRNA,Y_RNA,srpRNA,tRNA,snRNA,snoRNA,lncRNA,mRNA,tucpRNA \
    -o snakemake/mapping_small/sequential_mapping.snakemake
```

```bash
exseek.py --step mapping --dataset ${dataset}
```


## Quality control, adaptor removal and trimming

```bash
${exseek_path}/bin/exseek.py --step quality_control --dataset ${dataset}
```

Description of output files: [output_files](docs/output_files.md)

## Mapping (long RNA-seq)

```bash
${exseek_path}/bin/exseek.py --step mapping --dataset ${dataset}
```

Description of output files: [output_files](docs/output_files.md)


## Generate count matrix
```bash
${exseek_path}/bin/exseek.py --step count_matrix --dataset ${dataset}
```


**Count matrix**

* File path: `${output_dir}/count_matrix/transcript.txt`
* First row: sample IDs
* First column: feature names
* Feature name: `gene_id|gene_type|gene_name`


## Call domains for long RNA

```bash
${exseek_path}/bin/exseek.py --step call_domains --dataset ${dataset}
```


**Read count matrix**

* File path: `${output_dir}/count_matrix/domain_long.txt`
* First row: sample IDs
* First column: feature names
* Feature name: `gene_id|gene_type|gene_name|domain_id|transcript_id|start|end`

## Normalization

```bash
${exseek_path}/bin/exseek.py --step normalization --dataset ${dataset}
```
