## Call domains (from BAM files in genomic coordinates)
```bash
snakemake --snakefile snakemake/call_domains/Snakefile --configfile snakemake/config.json \
    --config 'output_dir'='output/scirep' \
    'sample_id_file'='metadata/sample_ids.scirep.txt'
```