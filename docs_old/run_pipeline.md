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

## Gini index comparison between exRNA and TCGA
```bash

for dataset in 'scirep' 'tcga_crc' 'tcga_hcc' 'GSE114711' 'GSE45722';do
    mkdir -p "output/${dataset}/analysis"
    echo bin/extract_bigwig.py abundant_rna_coverage_matrix \
        --matrix output/${dataset}/count_matrix/transcript.txt \
        --bigwig-pattern "output/${dataset}/tbigwig/{sample_id}.{gene_type}.bigWig" -n 100 \
        -o output/${dataset}/analysis/abundant_rna_coverage_matrix.h5
done
```

## Integrated dataset SciRep2016, GSE45722, GSE114711

### Copy files
```bash
dataset='exosome_small'
mkdir -p data/${dataset}
sources=(scirep GSE45722 GSE114711)
# sample_ids.txt
cat data/scirep/sample_ids.txt data/GSE45722/sample_ids.txt \
    data/GSE114711/sample_ids.txt > data/${dataset}/sample_ids.txt
# sample_classes.txt
{
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2}' data/scirep/sample_classes.txt 
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,"Healthy Control"}' data/GSE45722/sample_classes.txt
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,"Healthy Control"}' data/GSE114711/sample_classes.txt 
} > data/${dataset}/sample_classes.txt
# batch_info.txt
{
    echo -e 'sample_id\tpublication'
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,"scirep"}' data/scirep/sample_classes.txt 
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,"GSE45722_healthy"}' data/GSE45722/sample_classes.txt
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,"GSE114711_healthy"}' data/GSE114711/sample_classes.txt 
} > data/exosome_small/batch_info.txt

# count_matrix/domains_long.txt
mkdir -p output/${dataset}/count_matrix
bin/preprocess.py merge_data_frames \
    -i output/scirep/count_matrix/transcript.txt \
    -i output/GSE45722/count_matrix/transcript.txt \
    -i output/GSE114711/count_matrix/transcript.txt \
    --on feature --fillna 0 \
    -o output/${dataset}/count_matrix/transcript.txt
# domains_by_sample/20/05
mkdir -p output/${dataset}/domains_by_sample/20/05
rsync -ra output/scirep/domains_by_sample/20/05/*.bed output/${dataset}/domains_by_sample/20/05
rsync -ra output/GSE45722/domains_by_sample/20/05/*.bed output/${dataset}/domains_by_sample/20/05
rsync -ra output/GSE114711/domains_by_sample/20/05/*.bed output/${dataset}/domains_by_sample/20/05
# tbed
mkdir -p output/${dataset}/tbed
mkdir -p output/${dataset}/gbed

for d in ${sources[@]};do
    rsync -rav output/$d/tbed/ output/${dataset}/tbed
    rsync -rav output/$d/gbed/ output/${dataset}/gbed
done


```

### Analysis
```bash
bin/exseek.py call_domains -d exosome_small
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

## BigWig files
```bash
dest="/BioII/lulab_b/shared/projects/exSeek/output/exorbase/bigwig"
[ -d "$dest" ] || mkdir -p "$dest"
rsync -rav --delete /Share/home/shibinbin/projects/exSeek-dev/output/exorbase/bigwig/ "$dest/"
```

```bash
bin/create_igv.py -r genome/hg38 -g hg38 -i templates/igv/main.html \
     -o igv.html \
    --track output/exorbase/bigwig/SRR5679904.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5679904.genome_rmdup.+.bigWig \
    --track output/exorbase/bigwig/SRR5679905.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5679905.genome_rmdup.+.bigWig \
    --track output/exorbase/bigwig/SRR5679906.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5679906.genome_rmdup.+.bigWig \
    --track output/exorbase/bigwig/SRR5679907.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5679907.genome_rmdup.+.bigWig \
    --track output/exorbase/bigwig/SRR5679908.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5679908.genome_rmdup.+.bigWig \
    --track output/exorbase/bigwig/SRR5679909.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5679909.genome_rmdup.+.bigWig \
    --track output/exorbase/bigwig/SRR5687235.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5687235.genome_rmdup.+.bigWig \
    --track output/exorbase/bigwig/SRR5687236.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5687236.genome_rmdup.+.bigWig \
    --track output/exorbase/bigwig/SRR5687237.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5687237.genome_rmdup.+.bigWig \
    --track output/exorbase/bigwig/SRR5687238.genome_rmdup.-.bigWig \
    --track output/exorbase/bigwig/SRR5687238.genome_rmdup.+.bigWig
```

## Visualize domains
```bash
bin/report.py visualize_domains --sample-classes data/scirep/sample_classes.txt \
    --features output/scirep/feature_selection/filter.viper_count.Norm_CPM.Batch_RUV.domains_combined/Normal-CRC/random_forest.10.robust/features.txt \
    --genome-dir genome/hg38 \
    --flanking 20 -o tmp/visualize_domains.pdf --output-dir output/scirep
```

```bash
[ -d "output/scirep/visualize_domains" ] || mkdir -p output/scirep/visualize_domains
for features in $(cat /home/chenxupeng/projects/exseek/output/selected_feature/scirep/selected_features.txt);do
    compare_group=$(echo $features | cut -d'/' -f5 )
    bin/report.py visualize_domains --sample-classes data/scirep/sample_classes.txt \
        --features $features \
        --genome-dir genome/hg38 \
        --flanking 20 -o output/scirep/visualize_domains/${compare_group}.pdf --output-dir output/scirep
done
```


## Feature selection
```bash
bin/machine_learning.py  cross_validation \
    --matrix output/scirep/matrix_processing/filter.viper_count.Norm_CPM.Batch_RUV.domains_combined.txt \
    --sample-classes data/scirep/sample_classes.txt \
    -o tmp/cross_validation \
    --transpose \
    --positive-class 'Colorectal Cancer' --negative-class 'Healthy Control' \
    --cv-params '{"splitter": "stratified_shuffle_split", "n_splits": 5, "test_size": 0.2}' \
    --zero-fraction-filter \
    --log-transform '{"base": 2}' \
    --rpkm-filter '{"threshold": 10}' \
    --scaler robust \
    --selector robust --selector-params '{"cv": {"splitter": "stratified_shuffle_split", "n_splits": 10, "test_size": 0.2}}' \
    --grid-search \
    --grid-search-param-grid '{"n_estimators": [25, 50, 75], "max_depth": [3, 4, 5]}' \
    --grid-search-cv-params '{"splitter": "stratified_shuffle_split", "n_splits": 5, "test_size": 0.2}' \
    --classifier random_forest
```

## Notes

### 2018.12.29

* Feature selection with pairs of features
* Find enriched rRNA fragments
* Validate tRNA fragment features