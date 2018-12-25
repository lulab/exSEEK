## Generate sequential mapping snakefile
```bash
snakemake --snakefile snakemake/prepare_genome.snakemake --configfile snakemake/config.yaml --rerun-incomplete -k
```

## SciRep
```bash

bin/generate_snakemake.py sequential_mapping --rna-types rRNA,miRNA,piRNA,Y_RNA,srpRNA,tRNA,snRNA,snoRNA,lncRNA,mRNA,tucpRNA \
    -o snakemake/sequential_mapping.snakemake
snakemake --snakefile snakemake/mapping_small.snakemake --configfile config/scirep.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"


snakemake --snakefile snakemake/expression_matrix.snakemake --configfile config/scirep.yaml --rerun-incomplete -k

snakemake --snakefile snakemake/feature_selection.snakemake --configfile config/scirep.yaml --rerun-incomplete -k

snakemake --snakefile snakemake/feature_selection.snakemake --configfile config/scirep.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"

snakemake --snakefile snakemake/call_domains_long.snakemake --configfile config/scirep.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"
snakemake --snakefile snakemake/bigwig.snakemake --configfile config/scirep.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"

export PATH=$PWD/singularity/wrappers:$PATH
PATH="$PWD/singularity/wrappers:$PATH" snakemake --snakefile snakemake/normalization.snakemake --configfile config/scirep.yaml --rerun-incomplete -k

bin/report.py visualize_domains --sample-ids-file data/scirep/sample_ids.txt \
    --output-dir output/scirep \
    --count-matrix output/scirep/count_matrix/domains_combined.txt \
    --features output/scirep/feature_selection/filter.scimpute_count.Norm_CPM.Batch_RUV.domains_combined/Normal-CRC/random_forest.10.robust/features.txt \
    --chrom-sizes genome/hg38/chrom_sizes/transcriptome_genome \
    --output-file tmp/visualize_domains.pdf
bin/feature_selection.py calculate_clustering_score \
    --matrix output/scirep/matrix_processing/filter.scimpute_count.Norm_null.domains_combined.txt \
    --sample-classes data/scirep/sample_classes.txt --transpose

snakemake --snakefile snakemake/evaluate_features.snakemake --configfile config/scirep.yaml
```

## Lulab HCC
```bash
/Share/home/caojingyi/exRNA/process/18.new_hcc_lulab/Snakefile
snakemake --snakefile snakemake/mapping_small.snakemake --configfile config/lulab_hcc.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" \

snakemake --snakefile snakemake/mapping_small.snakemake --configfile config/lulab_hcc.yaml --rerun-incomplete -k
snakemake --snakefile snakemake/expression_matrix.snakemake --configfile config/lulab_hcc.yaml --rerun-incomplete -k
snakemake --snakefile snakemake/normalization.snakemake --configfile config/lulab_hcc.yaml --rerun-incomplete -k
snakemake --snakefile snakemake/feature_selection.snakemake --configfile config/lulab_hcc.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"

snakemake --snakefile snakemake/evaluate_features.snakemake --configfile config/lulab_hcc.yaml
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

snakemake --snakefile snakemake/mapping_long.snakemake --rerun-incomplete -k --configfile config/exorbase.yaml \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" -j60
snakemake --snakefile snakemake/normalization.snakemake --configfile config/exorbase.yaml --rerun-incomplete -k
snakemake --snakefile snakemake/feature_selection.snakemake --configfile config/exorbase.yaml --rerun-incomplete -k

```

## TCGA miRNA-seq (CRC)
```bash
snakemake --snakefile snakemake/bam_to_fastx.snakemake --rerun-incomplete -k --configfile config/tcga_crc.yaml 
snakemake --snakefile snakemake/mapping_small.snakemake --rerun-incomplete -k --configfile config/tcga_crc.yaml \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"
```

## TCGA miRNA-seq (HCC)
```bash
snakemake --snakefile snakemake/bam_to_fastx.snakemake --rerun-incomplete -k --configfile config/tcga_hcc.yaml 
snakemake --snakefile snakemake/mapping_small.snakemake --rerun-incomplete -k --configfile config/tcga_hcc.yaml \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"
snakemake --snakefile snakemake/bigwig.snakemake --configfile config/tcga_hcc.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"
snakemake --snakefile snakemake/call_domains_long.snakemake --configfile config/tcga_hcc.yaml --rerun-incomplete -k \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"
```

## GSE113994 (Healthy, cfRNA, PNAS 2018)
```bash
snakemake --snakefile snakemake/quality_control.snakemake --rerun-incomplete -k --configfile config/GSE113994.yaml \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"
snakemake --snakefile snakemake/mapping_small.snakemake --rerun-incomplete -k --configfile config/GSE113994.yaml \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}"
```

## GSE45722 (Healthy, exosome, BMC Genomics 2013)
```bash
snakemake --snakefile snakemake/quality_control.snakemake --rerun-incomplete -k --configfile config/GSE45722.yaml \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" -j40
snakemake --snakefile snakemake/mapping_small.snakemake --rerun-incomplete -k --configfile config/GSE45722.yaml \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" -j40
```

## GSE114711 (Healthy, exosome, Scientific Reports 2018)
```bash
snakemake --snakefile snakemake/quality_control.snakemake --rerun-incomplete -k --configfile config/GSE114711.yaml \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" -j40
snakemake --snakefile snakemake/mapping_small.snakemake --rerun-incomplete -k --configfile config/GSE114711.yaml \
    --cluster-config snakemake/cluster.yaml \
    --cluster "bsub -q {cluster.queue} -J {cluster.name} -e {cluster.stderr} -o {cluster.stdout} -R {cluster.resources} -n {cluster.threads}" -j40
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
ln -f -s /BioII/zhuyumin/exLocator/exosome_SR2018_GSE114711/1.fastq/*.fastq data/GSE114711/fastq
```

## Spike-in

### Small RNA-seq
`/BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/NEB/03.1811_T4PNK/ExiSEQ-spikeIn`

```bash
cp /BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/NEB/03.1811_T4PNK/ExiSEQ-spikeIn/ExiSEQ-spikeIn.fa genome/hg38/fasta/spikein_small.fa
samtools faidx genome/hg38/fasta/spikein_small.fa
```

### Long RNA-seq
`/BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/pico-smart/exSeek/ERCC-spikeIn`

```bash
cp /BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/pico-smart/exSeek/ERCC-spikeIn/ERCC92.fa genome/hg38/fasta/spikein_long.fa
samtools faidx genome/hg38/fasta/spikein_long.fa
```