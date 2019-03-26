### Input data files

| File name | Description |
| ------ | ----------- |
| `${input_dir}/fastq/${sample_id}.fastq` | Read files (single-end sequencing) |
| `${input_dir}/fastq/${sample_id}_1.fastq`, `${input_dir}/fastq/${sample_id}_2.fastq` | Read files (paired-end sequencing) |
| `${input_dir}/sample_ids.txt` | A text file with one sample ID per line. |
| `${input_dir}/sample_classes.txt` | A tab-deliminated file (with header) with two columns: sample_id, label (optional) |
| `${input_dir}/batch_info.txt` | A comma-deliminated file (with header) with at least two columns: sample_id, batch1, batch2, ... (optional) |
| `${input_dir}/compare_groups.yaml` | A YAML file defining positive and negative classes. (optional) |
| `${config_dir}/${dataset}.yaml` | A YAML file for configuration parameters for the dataset |

Example config file: `config/long_pe_example.yaml`


## Quality control (before adapter removal)

```bash
exseek.py quality_control -d ${dataset}
```

## Remove adapter

```bash
exseek.py cutadapt -d ${dataset}
```

## Start with clean reads

```bash
mkdir -p ${output_dir}/cutadapt
ln -f -s ${fastq_dir}/*.fastq.gz ${output_dir}/cutadapt
```

## Quality control (after adapter removal)

```bash
exseek.py quality_control_clean -d ${dataset}
```

## Rename fastq files

```bash
exseek.py rename_fastq -d ${dataset}
```

## Mapping

```bash
exseek.py mapping -d ${dataset}
```

## Count matrix

```bash
exseek.py count_matrix -d ${dataset}
```
