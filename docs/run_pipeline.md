## Generate sequential mapping snakefile
```bash
snakemake --snakefile snakemake/prepare_genome.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k
```

## SciRep
```bash

bin/generate_snakemake.py sequential_mapping --rna-types rRNA,miRNA,piRNA,Y_RNA,srpRNA,tRNA,snRNA,snoRNA,lncRNA,mRNA,tucpRNA \
    -o snakemake/sequential_mapping.snakemake
snakemake --snakefile snakemake/mapping_small.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" output_dir="output/scirep" \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"

snakemake --snakefile snakemake/mapping_small.snakemake \
    --cluster-config snakemake/cluster.yaml \
    --configfile snakemake/config.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" sample_id_file="metadata/sample_ids.scirep.txt" \
    output_dir="output/scirep" \
    --rerun-incomplete -k -j60 

snakemake --snakefile snakemake/expression_matrix.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" output_dir="output/scirep"

snakemake --snakefile snakemake/feature_selection.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" output_dir="output/scirep" \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"

snakemake --snakefile snakemake/call_domains_long.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" output_dir="output/scirep" \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"


export PATH=$PWD/singularity/wrappers:$PATH
snakemake --snakefile snakemake/normalization.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/scirep" output_dir="output/scirep"
```

## Lulab HCC
```bash
/Share/home/caojingyi/exRNA/process/18.new_hcc_lulab/Snakefile
snakemake --snakefile snakemake/mapping_small.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc" output_dir="output/lulab_hcc" -j8

snakemake --snakefile snakemake/mapping_small.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k \
    --configfile snakemake/config.yaml \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc" output_dir="output/lulab_hcc" -j8

snakemake --snakefile snakemake/expression_matrix.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k \
    --configfile snakemake/config.yaml \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc" output_dir="output/lulab_hcc" -j8

snakemake --snakefile snakemake/call_domains_long.snakemake --rerun-incomplete -k \
    --configfile snakemake/config.yaml \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc"  output_dir="output/lulab_hcc" -j8

snakemake --snakefile snakemake/normalization.snakemake --rerun-incomplete -k --configfile snakemake/config.yaml \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc"  output_dir="output/lulab_hcc" -j8

snakemake --snakefile snakemake/feature_selection.snakemake --rerun-incomplete -k --configfile snakemake/config.yaml \
    --config adaptor="AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC" data_dir="data/lulab_hcc"  output_dir="output/lulab_hcc" -j8
```

## exoRBase
```bash
# subsample reads
[ -d "data/exorbase_test/fastq" ] || mkdir -p "data/exorbase_test/fastq"
cp data/exorbase/*.txt data/exorbase_test
for f in data/exorbase/fastq/*.fastq;do
    sample_id=$(basename $f)
    sample_id=${sample_id/.fastq/}
    echo data/exorbase_test/fastq/${sample_id}.fastq
    head -n 400000 $f > data/exorbase_test/fastq/${sample_id}.fastq
done
snakemake --snakefile snakemake/mapping_long.snakemake --rerun-incomplete -k --configfile snakemake/config.yaml \
    --config data_dir="data/exorbase_test"  output_dir="output/exorbase_test"


snakemake --snakefile snakemake/mapping_long.snakemake --rerun-incomplete -k --configfile snakemake/config.yaml \
    --config data_dir="data/exorbase"  output_dir="output/exorbase" \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" -j60

```

## Data

### exoRBase

```bash
[ -d data/exorbase/fastq ] || mkdir -p data/exorbase/fastq
ln -f -s /BioII/lulab_b/shared/projects/exRNA/published_exRNA/exosome_exoRBase/exosome_GSE100063_CRC/fastq/*.fastq data/exorbase/fastq
ln -f -s /BioII/lulab_b/caojingyi/exoRBase/fastq/HCC/*.fastq data/exorbase/fastq
ln -f -s /BioII/lulab_b/caojingyi/exoRBase/fastq/Normal/*.fastq  data/exorbase/fastq
ln -f -s /BioII/lulab_b/shared/projects/exRNA/published_exRNA/exosome_exoRBase/exosome_GSE100232_PAAD/fastq/*.fastq data/exorbase/fastq
ln -f -s /BioII/lulab_b/shared/projects/exRNA/published_exRNA/exosome_exoRBase/exosome_GSE99985_CHD/fastq/*.fastq data/exorbase/fastq
ln -f -s /BioII/lulab_b/caojingyi/exoRBase/fastq/*.fastq data/exorbase/fastq
ls data/exorbase/fastq | cut -d'_' -f1 | sort | uniq > data/exorbase/sample_ids.txt
```