# Small RNA-seq mapping

## Create spikein directory

```bash
# output spike-in directory
spikein_dir="data/${dataset}/spikein"
# input spike-in sequences
spikein_fa="data/${dataset}/spikein.fa"
# create sub-directories
mkdir -p "$spikein_dir"
for sub_dir in fasta bed chrom_sizes transcript_table index/bowtie2;do
    mkdir -p "$spikein_dir/$sub_dir"
done
# create files
cp "$spikein_fa" "$spikein_dir/fasta/spikein.fa"
samtools faidx "$spikein_dir/fasta/spikein.fa"
cut -f1,2 "$spikein_dir/fasta/spikein.fa.fai" > "$spikein_dir/chrom_sizes/spikein"
{
    echo -e 'chrom\tstart\tend\tname\tscore\tstrand\tgene_id\ttranscript_id\tgene_name\ttranscript_name\tgene_type\ttranscript_type\tsource'
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,0,$2,$1,0,"+",$1,$1,$1,$1,"spikein","spikein","spikein"}' "$spikein_dir/fasta/spikein.fa.fai"
} > "$spikein_dir/transcript_table/spikein.txt"
bowtie2-build "$spikein_dir/fasta/spikein.fa" "$spikein_dir/index/bowtie2/spikein"
```

## Update sequential mapping order

The default mapping order is set as `rna_type` variable in `snakemake/default_config.yaml`:

```yaml
rna_types: [rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, 
  snRNA, srpRNA, tRNA, tucpRNA, Y_RNA]
```

You can change the mapping order by add a `rna_type` variable in `config/${dataset}.yaml`. For example, add spike-in sequences as the first RNA type:

```yaml
rna_types: [spikein, rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, 
  snRNA, srpRNA, tRNA, tucpRNA, Y_RNA]
```

```bash
exseek.py update_sequential_mapping -d ${dataset}
```

## Add new reference sequence

If a new RNA type is added, you should also add a sequence file in FASTA format: `${genome_dir}/fasta/${rna_type}.fa`. Then build a FASTA index \(`${genome_dir}/fasta/${rna_type}.fa.fai`\):

```bash
samtools faidx ${genome_dir}/fasta/${rna_type}.fa
```

Then build a bowtie2 index \(`${genome_dir}/index/bowtie2/${rna_type}`\):

```bash
bowtie2-build ${genome_dir}/fasta/${rna_type}.fa ${genome_dir}/index/bowtie2/${rna_type}
```

## Quality control \(before adapter removal\)

```bash
exseek.py quality_control -d ${dataset}
```

## Remove adapter

```bash
exseek.py cutadapt -d ${dataset}
```

## Start clean reads

```bash
mkdir -p ${output_dir}/cutadapt
ln -f -s ${fastq_dir}/*.fastq.gz ${output_dir}/cutadapt
# convert to fasta
exseek.py fastq_to_fasta -d ${dataset}
```

## Quality control \(after adapter removal\)

```bash
exseek.py quality_control_clean -d ${dataset}
```

## Mapping

```bash
exseek.py mapping -d ${dataset}
```

## Generate BigWig files

```bash
exseek.py bigwig -d ${dataset}
```

## Call domains

```bash
exseek.py call_domains -d ${dataset}
```

## Count matrix

```bash
exseek.py count_matrix -d ${dataset}
```

## Combine domains with small RNA

```bash
exseek.py combine_domains -d ${dataset}
```

