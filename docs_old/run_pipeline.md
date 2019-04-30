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

## GSE112343 (Healthy, cfRNA, bioRxiv 2018)
```bash
awk 'BEGIN{OFS="\t";FS="\t";print "sample_id\tlabel"}{print $0,"Healthy"}' data/GSE112343/sample_ids.txt > data/GSE112343/sample_classes.txt
```

## GSE113994 (Healthy, cfRNA, PNAS 2018)

### Get adapters:
```bash
curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE113994&targ=self&form=text&view=quick" \
    | awk '{match($0, /!Series_sample_id = (GSM[0-9]+)/, a); if(a[1])print a[1]}' > data/GSE113994/gsm_ids.txt
{
for gsm_id in $(cat data/GSE113994/gsm_ids.txt);do
    adapter=$(curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm_id}&targ=self&form=text&view=quick" \
        | grep '^!Sample_characteristics_ch1' | grep adapter \
        | awk '{match($0,/\(([ATCG]+)\)/,a);print a[1]}')
    echo -e "$gsm_id\t$adapter"
done
} > data/GSE113994/adapters_by_gsm_id.txt
awk 'BEGIN{OFS="\t";FS="\t"} (FNR==NR){if(NR>1){srr[$7]=$5};next} {print srr[$1],$2}' \
    data/GSE113994/SraRunTable.txt data/GSE113994/adapters_by_gsm_id.txt \
    > data/GSE113994/adapters.txt
```

## GSE53080 (Healthy, cfRNA, PNAS 2014)
```bash
curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE53080&targ=self&form=text&view=quick" \
    | awk '{match($0, /!Series_sample_id = (GSM[0-9]+)/, a); if(a[1])print a[1]}' > data/GSE53080/gsm_ids.txt
{
for gsm_id in $(cat data/GSE53080/gsm_ids.txt);do
    adapter=$(curl -s "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=${gsm_id}&targ=self&form=text&view=quick" \
        | grep '^!Sample_description' \
        | awk "{match(\$0,/complete 3'-adapter: ([ATCG]+)/,a);if(a[1]) print a[1]}")
    echo -e "$gsm_id\t$adapter"
done
} > data/GSE53080/adapters_by_gsm_id.txt
awk 'BEGIN{OFS="\t";FS="\t"} (FNR==NR){if(NR>1){srr[$11]=$9};next} {print srr[$1],$2}' \
    data/GSE53080/SraRunTable.txt data/GSE53080/adapters_by_gsm_id.txt \
    > data/GSE53080/adapters.txt
awk -F'\t' '$24=="Plasma"&&$14=="None"{print $9}' data/GSE53080/SraRunTable.txt > data/GSE53080/sample_ids.txt
awk 'BEGIN{OFS="\t";FS="\t";print "sample_id\tlabel"}{print $0,"Healthy"}' \
    data/GSE53080/sample_ids.txt > data/GSE53080/sample_classes.txt
# Spike-ins: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3508525/
data/GSE53080/spikein/
```

## GSE67489 (healthy, cfRNA, Genome Medicine)
```bash

```


## Combine cfRNA datasets

### Prepare input files
```bash
datasets="lulab_hcc GSE94582 GSE53080 GSE113994"
healthy_datasets="GSE94582 GSE53080 GSE113994"
pvalue=05
bin_size=20

# directories
for d in "data/cfRNA" "output/cfRNA";do
    [ -d "$d" ] || mkdir -p "$d"
done
# sample_ids
{
for dataset in $datasets;do
    cat data/$dataset/sample_ids.txt
done
} > data/cfRNA/sample_ids.txt
# sample_classes
{
cat data/lulab_hcc/sample_classes.txt
for dataset in $healthy_datasets;do
    awk 'BEGIN{OFS="\t"}{print $1,"Normal"}' data/$dataset/sample_ids.txt
done 
} > data/cfRNA/sample_classes.txt
# batch info
{
echo -e "sample_id\tdataset"
for dataset in $datasets;do
    awk -v d=$dataset 'BEGIN{OFS="\t"}{print $1,d}' data/$dataset/sample_ids.txt
done
} > data/cfRNA/batch_info.txt

for d in "counts/transcript" "domains_by_sample/$bin_size/$pvalue" "tbed_long_RNA";do
    [ -d "output/cfRNA/$d" ] || mkdir -p "output/cfRNA/$d"
done

# link counts
for dataset in $datasets;do
    for sample_id in $(cat data/$dataset/sample_ids.txt);do
        # transcript counts
        ln -f -r -s output/$dataset/counts/transcript/$sample_id output/cfRNA/counts/transcript/$sample_id
        # domains_by_sample
        ln -f -r -s output/$dataset/domains_by_sample/$bin_size/$pvalue/${sample_id}.bed output/cfRNA/domains_by_sample/$bin_size/$pvalue/${sample_id}.bed
        # tbed_long_RNA
        ln -f -r -s output/$dataset/tbed_long_RNA/${sample_id}.bed.gz output/cfRNA/tbed_long_RNA/${sample_id}.bed.gz
    done
done

# link BAM files and BigWig files
for dataset in $datasets;do
    for d in "gbam_sorted" "gbam" "bigwig" "bigwig_normalized" "tbigwig" "tbigwig_normalized" "bigwig_normalized_log2";do
        [ -d "output/cfRNA/$d" ] || mkdir -p "output/cfRNA/$d"
        ln -f -r -s output/$dataset/$d/* output/cfRNA/$d/
    done
done
```

## GSE94582 (healthy, cfRNA, Nature Biotechnology)
```bash
awk 'BEGIN{OFS="\t"}{if($1 ~ /^4N/){print $1,"NNNNTGGAATTCTCGGGTGCCAAGG"}else{print $1,"TGGAATTCTCGGGTGCCAAGG"}}' \
    data/GSE94582/sample_ids.txt > data/GSE94582/adapters.txt

awk 'BEGIN{OFS="\t"}{if($1 ~ /^4N/){print $1,"GTTCAGAGTTCTACAGTCCGACGATCNNNN"}else{print $1,"GTTCAGAGTTCTACAGTCCGACGATC"}}' \
    data/GSE94582/sample_ids.txt > data/GSE94582/adapters_5p.txt

awk 'BEGIN{OFS="\t";FS="\t";print "sample_id\tlabel"}{print $0,"Healthy"}' \
    data/GSE94582/sample_ids.txt > data/GSE94582/sample_classes.txt

data_dir=/BioII/lulab_b/shared/projects/exRNA/published_exRNA/cfRNA_NatBio2018_GSE94582_Healthy/01.preprocess
mkdir -p output/GSE94582/cutadapt
for sample_id in A_Lab5 B_Lab2 B_Lab4 B_Lab5 B_Lab6 C_Lab6 D_Lab1 Xu_Lab5;do
    pigz -c ${data_dir}/N4_${sample_id}.cutadapt_2.fastq > output/GSE94582/cutadapt/4N_${sample_id}.fastq.gz
done
# remove barcode
for sample_id in A_Lab5 B_Lab2 B_Lab4 B_Lab5 B_Lab6 C_Lab6 D_Lab1 Xu_Lab5;do
    cutadapt -u 4 -u -4 -m 16 --trim-n -q 30 ${data_dir}/N4_${sample_id}.cutadapt_1.fastq \
        | pigz -c -p 4 > output/GSE94582/cutadapt/4N_${sample_id}.fastq.gz
done
for sample_id in CleanTag_Lab5 NEBNext_Lab1 NEBNext_Lab2 NEBNext_Lab3 NEBNext_Lab4 NEBNext_Lab5 TruSeq_Lab1 TruSeq_Lab2 TruSeq_Lab3 TruSeq_Lab4 TruSeq_Lab5 TruSeq_Lab6;do
    pigz -c ${data_dir}/${sample_id}.cutadapt.fastq > output/GSE94582/cutadapt/${sample_id}.fastq.gz
done
```

### Create spike-in directory
```bash
spikein_dir="data/GSE113994/spikein"
spikein_fa="tmp/spikein.fa"

mkdir -p "$spikein_dir"
mkdir -p "$spikein_dir/fasta"
cp "$spikein_fa" "$spikein_dir/fasta/spikein.fa"
samtools faidx "$spikein_dir/fasta/spikein.fa"
mkdir -p "$spikein_dir/bed"
mkdir -p "$spikein_dir/chrom_sizes"
mkdir -p "$spikein_dir/transcript_table"
cut -f1,2 "$spikein_dir/fasta/spikein.fa.fai" > "$spikein_dir/chrom_sizes/spikein"
{
    echo -e 'chrom\tstart\tend\tname\tscore\tstrand\tgene_id\ttranscript_id\tgene_name\ttranscript_name\tgene_type\ttranscript_type\tsource'
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,0,$2,$1,0,"+",$1,$1,$1,$1,"spikein","spikein","spikein"}' "$spikein_dir/fasta/spikein.fa.fai"
} > "$spikein_dir/transcript_table/spikein.txt"
mkdir -p "$spikein_dir/index/bowtie2"
bowtie2-build "$spikein_dir/fasta/spikein.fa" "$spikein_dir/index/bowtie2/spikein"
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
awk 
# sample_classes.txt
{
    echo -e 'sample_id\tlabel'
    awk 'BEGIN{OFS="\t";FS="\t"}NR>1{print $1,$2}' data/scirep/sample_classes.txt 
    awk 'BEGIN{OFS="\t";FS="\t"}NR>1{print $1,"Healthy Control"}' data/GSE45722/sample_classes.txt
    awk 'BEGIN{OFS="\t";FS="\t"}NR>1{print $1,"Healthy Control"}' data/GSE114711/sample_classes.txt 
} > data/${dataset}/sample_classes.txt
# batch_info.txt
{
    echo -e 'publication'
    awk 'BEGIN{OFS="\t";FS="\t"}NR>1{print $1,"scirep"}' data/scirep/sample_classes.txt 
    awk 'BEGIN{OFS="\t";FS="\t"}NR>1{print $1,"GSE45722_healthy"}' data/GSE45722/sample_classes.txt
    awk 'BEGIN{OFS="\t";FS="\t"}NR>1{print $1,"GSE114711_healthy"}' data/GSE114711/sample_classes.txt 
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

```bash
Rscript /Share/home/shibinbin/projects/exSeek-dev/bin/matrix-process.R -s batch_removal -i output/exosome_small/count_matrix/domains_combined.txt -c data/exosome_small/sample_classes.txt -b data/exosome_small/batch_info.txt --filterout output/exosome_small/matrix_processing/ --imputeout output/exosome_small/matrix_processing/ --normalizeout output/exosome_small/matrix_processing/ --batchremoveout output/exosome_small/matrix_processing/ --imputemethod null --filtercount 5 --filtersample 10 --imputecluster 5 -p 1 --normmethod TMM --normtopk 20 --removetype miRNA,piRNA --cvthreshold 0.5 --batchmethod Combat --batchindex 1
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
    --cv-params '{"splitter": "stratified_shuffle_split", "n_splits": 1, "test_size": 0.2}' \
    --zero-fraction-filter --zero-fraction-filter-threshold 0.8 \
    --log-transform --log-transform-base 2 \
    --rpm-filter --rpm-filter-threshold 10 \
    --fold-change-filter --fold-change-filter-direction 'up' \
    --scaler robust --scaler-params '{"with_centering": false}' \
    --selector robust --selector-params '{"cv": {"splitter": "stratified_shuffle_split", "n_splits": 10, "test_size": 0.2}}' \
    --grid-search \
     --grid-search-params '{"cv": {"splitter": "stratified_shuffle_split", "n_splits": 5, "test_size": 0.2},
     "param_grid": {"C": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2, 1e3, 1e4, 1e5]}}' \
    --classifier logistic_regression

    --grid-search-params '{"cv": {"splitter": "stratified_shuffle_split", "n_splits": 5, "test_size": 0.2},
     "param_grid": {"C": [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 1e2, 1e3, 1e4, 1e5]}}' \
    --classifier linear_svm

    --grid-search-params '{"cv": {"splitter": "stratified_shuffle_split", "n_splits": 5, "test_size": 0.2},
     "param_grid": {"n_estimators": [25, 50, 75], "max_depth": [3, 4, 5]}}' \
    --classifier random_forest
```

## Notes

### 2018.12.29

* Feature selection with pairs of features
* Find enriched rRNA fragments
* Validate tRNA fragment features

## Pico 3V3

### Prepare input clean reads
```bash
/BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/pico-smart/05.pico-3v3/01.cutadapt/*.filter.fastq
[ -d "output/pico_3v3/cutadapt" ] || mkdir -p "output/pico_3v3/cutadapt"
for sample_id in $(cat data/pico_3v3/sample_ids.txt);do
    for mate in 1 2;do
        echo "gzip ${sample_id}_${mate}"
        pigz -c /BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/pico-smart/05.pico-3v3/01.cutadapt/${sample_id}_${mate}.filter.fastq \
            > "output/pico_3v3/cutadapt/${sample_id}_${mate}.fastq.gz"
    done
done
```

### Run pipeline
```bash
bin/exseek.py quality_control 
```

## Known biomarkers
```bash
for dataset in "scirep" "exorbase";do
    for compare_group in $(ls output/$dataset/evaluate_features/matrix);do
        [ -d "data/$dataset/known_biomarkers/$compare_group" ] || mkdir -p "data/$dataset/known_biomarkers/$compare_group"
        for matrix_file in $(ls output/$dataset/evaluate_features/matrix/$compare_group);do
            awk 'NR>1{print $1}' output/$dataset/evaluate_features/matrix/$compare_group/$matrix_file \
                > data/$dataset/known_biomarkers/$compare_group/$matrix_file
        done
    done
done
for compare_group in "Normal-CRC_S1" "Normal-CRC_S2" "Normal-CRC_S3" "Normal-CRC_S4";do
    cp -r "data/scirep/known_biomarkers/Normal-CRC" "data/scirep/known_biomarkers/$compare_group"
done
```

## Differential expression
```bash
for method in edger_exact edger_glmlrt edger_glmqlf wilcox;do
    echo "Method: $method"
    bin/differential_expression.R -i output/scirep/matrix_processing/filter.domains_combined.txt \
        --method "$method" \
        -c data/scirep/sample_classes.txt \
        -p 'Colorectal Cancer Stage 1' \
        -n 'Healthy Control' \
        -o tmp/${method}.txt
done

```

## Quake 2018

### Prepare input files
```bash
mkdir -p data/quake_2018
mkdir -p output/quake_2018/cutadapt
for sample_id in $(cat data/quake_2018/sample_ids.txt);do
    for mate in 1 2;do
        ln -f -s /BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/pico-smart/Quake-2018-science.exSeek/01.preprocess/${sample_id}_${mate}.filter.fastq \
        output/quake_2018/cutadapt/${sample_id}_${mate}.fastq
    done
done
```

## CRC NEB 

### Prepare clean reads
```bash
source_dir="/BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/NEB/cutadapt"
pigz -c $source_dir/Colon1-qiazol_combined_R1.cutadapt.fastq > Colon1-qiazol_combined.fastq.gz
pigz -c $source_dir/Colon2-qiazol_combined_R1.cutadapt.fastq > Colon2-qiazol_combined.fastq.gz
pigz -c $source_dir/crc-2375612.cutadapt.fastq > crc-2375612.fastq.gz
pigz -c $source_dir/crc-4.cutadapt.fastq > crc-4.fastq.gz
cp $source_dir/crc-2382274.cutadapt.fastq.gz crc-2382274.fastq.gz
cp $source_dir/crc-2382274_T4PNK.cutadapt.fastq.gz crc-2382274_T4PNK.fastq.gz
cp $source_dir/crc-2382711.cutadapt.fastq.gz crc-2382711.fastq.gz
cp $source_dir/crc-2382711_T4PNK.cutadapt.fastq.gz crc-2382711_T4PNK.fastq.gz
```

## CRC S-smart

### Prepare clean reads
```bash
source_dir="/BioII/lulab_b/wangsiqi/exRNA/exRNA-panel/s-smart/cutadapt"
pigz -c $source_dir/SMART-sm11_1.cutadapt.fastq > SMART-sm11_1.fastq.gz
pigz -c $source_dir/SMART-sm9_1.cutadapt.fastq > SMART-sm9_1.fastq.gz
pigz -c $source_dir/smart_1.cutadapt_2.fastq > smart_1.fastq.gz
pigz -c $source_dir/T4PNK-smart_1.cutadapt_2.fastq > T4PNK-smart_1.fastq.gz
```

## RNA Structure comparison

### get sequences
```bash

```

```
/Share/home/shibinbin/projects/exSeek-dev/bin/machine_learning.py cross_validation --matrix output/cfRNA/matrix_processing/filter.null.Norm_TMM.Batch_limma_1.mirna_and_domains_rna.txt --sample-classes data/cfRNA/sample_classes.txt --output-dir output/cfRNA/cross_validation/filter.null.Norm_TMM.Batch_limma_1.mirna_and_domains_rna/Normal-stage_A/logistic_regression.6.max_features.any --transpose --positive-class stage_A --negative-class Normal --cv-params '{"splitter": "stratified_shuffle_split", "n_splits": 50, "test_size": 0.1}' --log-transform --log-transform-params '{"base": 2, "pseudo_count": 4}' --scaler zscore --scaler-params '{"with_mean": true}' --selector max_features --selector-params '{}' --n-features-to-select 6 --grid-search --grid-search-params '{"param_grid": {"logistic_regression": {"C": [1e-05, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]}, "logistic_regression_l1": {"C": [1e-05, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]}, "logistic_regression_l2": {"C": [1e-05, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]}, "linear_svm": {"C": [1e-05, 0.0001, 0.001, 0.01, 0.1, 1, 10, 100, 1000, 10000, 100000]}, "random_forest": {"n_estimators": [25, 50, 75], "max_depth": [3, 4, 5]}, "decision_tree": {"max_depth": [2, 3, 4, 5, 6, 7, 8]}}, "cv": {"splitter": "stratified_kfold", "n_splits": 5}, "iid": false, "scoring": "roc_auc"}' --sample-weight auto --classifier logistic_regression --classifier-params '{"penalty": "l2", "solver": "liblinear"}'
```