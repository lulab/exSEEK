## Generate sequential mapping snakefile
```bash
snakemake --snakefile snakemake/prepare_genome/Snakefile --configfile snakemake/config.json --rerun-incomplete -k
snakemake --snakefile snakemake/prepare_genome/Snakefile \
    --cluster-config snakemake/cluster.json \
    --configfile snakemake/config.json \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" \
    --rerun-incomplete -k -j50
```

```bash

bin/generate_snakemake.py sequential_mapping --rna-types rRNA,miRNA,piRNA,Y_RNA,srpRNA,tRNA,snRNA,snoRNA,lncRNA,mRNA,tucpRNA \
    -o snakemake/mapping_small/sequential_mapping
snakemake --snakefile snakemake/mapping_small/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k

snakemake --snakefile snakemake/expression_matrix/Snakefile --configfile snakemake/config.yaml --rerun-incomplete -k
```