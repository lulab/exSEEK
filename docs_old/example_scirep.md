## Raw data
```bash
[ -d "data/scirep/fastq" ] || mkdir -p "data/scirep/fastq"
for sample_id in $(cat metadata/sample_ids.scirep.txt);do
    ln -s -f /BioII/lulab_b/shared/projects/exRNA/published_exRNA/exosome_SR2017_GSE71008/00.rawdata/${sample_id}.fastq data/scirep/fastq/${sample_id}.fastq
done
```

## Call domains (from BAM files in genomic coordinates)
```bash
snakemake --snakefile snakemake/call_domains/Snakefile --configfile snakemake/config.json \
    --config 'output_dir'='output/scirep' \
    'sample_id_file'='metadata/sample_ids.scirep.txt'
```