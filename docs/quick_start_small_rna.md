## Installation

```bash
git clone https://github.com/ltbyshi/exseek-pipeline.git
cd exseek-pipeline/
python setup.py develop
```

## Small RNA

* **Input directory**: `data/example_small`
* **Configuration file**: `config/example_small.yaml`
* **Genome files**: `genome/hg38/fasta`, `genome/hg38/gtf_longest_transcript`, `genome/hg38/bed`

```bash
exseek -d example_small build_index
exseek -d example_small quality_control
exseek -d example_small cutadapt
exseek -d example_small quality_control_clean
exseek -d example_small mapping
exseek -d example_small bigwig
exseek -d example_small call_domains
exseek -d example_small count_matrix
exseek -d example_small normalization
exseek -d example_small feature_selection
```

