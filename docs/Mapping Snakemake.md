# Mapping Snakemake
## Check Validity of Sample_ids
check whether sample IDs contain '.' or not
## Check Validity of RNA_types
check whether RNA types contain '.' or not
## Parameter Settings
| ds| ds |
| :--- | :--- |
|`data_dir`|path to `.fastq` files|
|`data_dir`|path to `.fastq` files|
|`rna_types`|map to different types of RNA indexes|
|`adaptor`|reads have adaptors, and software cut them off provided with sequences |
|`min_read_length`| filter too short reads |
|`genome_dir`|  |
|`max_read_length`| filter too long reads |
|`min_base_quality`| base quality control |
|`temp_dir`| store temporary files |

## mapping statistics
###read\_counts_raw
Count reads in `.fastq` files of raw data

```
wc -l < {input} | awk '{{print int($0/4)}}' > {output}
```
###read\_counts_mapped
Count reads in `.bam` files of mapped reads

```
bamtools count -in {input} > {output}
```
###read\_counts_unmapped
Count reads in `.fa.gz` files of mapped reads

```
pigz -p {threads} -d -c {input} | wc -l | awk '{{print int($0/2)}}' > {output}
```
###summarize\_read_counts

###mapped\_read_length
run python script to count reads length of different `.bam` files as outputs of sequential mapping

```
bin/statistics.py read_length_hist --max-length 600 -i {input} -o {output}
```
###merge\_mapped\_read_length


###fastqc
fastqc of raw data

### parse\_fastqc_data
### summarize\_fastqc_ipynb
### summarize\_fastqc_html

### cutadapt
`cutadapt`: cutadapt removes adapter sequences from high-throughput sequencing reads.
```
cutadapt -a {params.adaptor} -m {params.min_read_length} --trim-n -q {min_base_quality}          --too-short-output >(pigz -c -p {threads} > {output.too_short}) -o {output.trimmed} {input}
```
### fastq\_to_fasta
Change file attributes to remove quality information

### tbam\_to_gbam
convert transcript coordinate BAM alignments file into a genomic coordinate BAM alignments file

```
rsem-tbam2gbam {params.index} {input.bam} 
```

### sort_gbam
```
samtools sort {input} > {output.bam}
samtools index {output.bam}
```
### gbam\_to_bedgraph

### gbedgraph\_to_bigwig

### sort_tbam
```
samtools sort -T {params.temp_dir} -o {output} {input}
```

### collect\_alignment\_summary_metrics
Produces a summary of alignment metrics from a SAM or BAM file.

```
picard CollectAlignmentSummaryMetrics I={input} O={output}
```
### count\_reads_intron
Provided with `.bed` file containing intron loci and `other.bam` file containing reads mapped to hg38, report overlaps to retrieve intron stats.

```
bedtools intersect -wa -s -a {input.bam} -b {input.bed} | wc -l > {output}
```

### count\_reads_promoter
Provided with `.bed` file containing promoter loci and `other.bam` file containing reads mapped to hg38, report overlaps to retrieve promoter stats.

### count\_reads_enhancer
Provided with `.bed` file containing enhancer loci and `other.bam` file containing reads mapped to hg38, report overlaps to retrieve enhancer stats.

### map_circRNA
The software aligns unmapped reads to cicrRNA index
...
```
pigz -d -c other.fa.gz | bowtie2 -f -p {threads} --norc --sensitive --no-unal --un-gz circRNA.aligner.fa.gz -x circRNA - -S - | bin/preprocess.py filter_circrna_reads --filtered-file >(samtools view -b -o {output.bam_filtered}) | samtools view -b -o {output.bam}
```