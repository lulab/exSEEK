## Generate sequential mapping snakefile
```bash
snakemake --snakefile snakemake/prepare_genome/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k
snakemake --snakefile snakemake/prepare_genome/Snakefile \
    --cluster-config snakemake/cluster.json \
    --configfile snakemake/config.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" \
    --rerun-incomplete -k -j40
```

## SciRep
```bash

bin/generate_snakemake.py sequential_mapping --rna-types rRNA,miRNA,piRNA,Y_RNA,srpRNA,tRNA,snRNA,snoRNA,lncRNA,mRNA,tucpRNA \
    -o snakemake/mapping_small/sequential_mapping
snakemake --snakefile snakemake/mapping_small/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" sample_id_file="metadata/sample_ids.scirep.txt" \
    output_dir="output/scirep" 
snakemake --snakefile snakemake/mapping_small/Snakefile \
    --cluster-config snakemake/cluster.json \
    --configfile snakemake/config.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" sample_id_file="metadata/sample_ids.scirep.txt" \
    output_dir="output/scirep" \
    --rerun-incomplete -k -j60 

snakemake --snakefile snakemake/expression_matrix/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" sample_id_file="metadata/sample_ids.scirep.txt" \
    output_dir="output/scirep"

snakemake --snakefile snakemake/call_domains_long/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" sample_id_file="metadata/sample_ids.scirep.txt" \
    output_dir="output/scirep" \
    --cluster-config snakemake/cluster.json \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" \
    
```

## Lulab HCC
```bash
/Share/home/caojingyi/exRNA/process/18.new_hcc_lulab/Snakefile
snakemake --snakefile snakemake/mapping_small/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.json \
    --configfile snakemake/config.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc" sample_id_file="metadata/sample_ids.lulab_hcc.txt" \
    output_dir="output/lulab_hcc" -j8
snakemake --snakefile snakemake/mapping_small/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k \
    --configfile snakemake/config.yaml \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc" sample_id_file="metadata/sample_ids.lulab_hcc.txt" \
    output_dir="output/lulab_hcc" -j8

snakemake --snakefile snakemake/expression_matrix/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k \
    --configfile snakemake/config.yaml \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc" sample_id_file="metadata/sample_ids.lulab_hcc.txt" \
    output_dir="output/lulab_hcc" -j8

snakemake --snakefile snakemake/call_domains_long/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k \
    --configfile snakemake/config.yaml \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc" sample_id_file="metadata/sample_ids.lulab_hcc.txt" \
    output_dir="output/lulab_hcc" -j8

```