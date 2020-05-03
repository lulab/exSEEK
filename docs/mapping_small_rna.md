# Mapping small RNA-seq

## Prepare genome annotation

For mapping of small RNA-seq reads, exSEEK adopts sequential mapping strategy, which assign reads to gene annotations sequentially according the the ordered defined by the user.
By default, exSEEK assign reads in the following order:

spike-in, rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, snRNA, srpRNA, tRNA, tucpRNA, Y_RNA, genome, circRNA

We derived the genome annotation file from various sources: 

| Type | Number of genes | Source |
| :--- | :--- | :--- |
| miRNA | 1917 | miRBase hairpin \(Version 22\) |
| piRNA | 23431 | piRNABank |
| lncRNA | 15778 | GENCODE V27 and mitranscriptome |
| rRNA | 37 | NCBI refSeq 109 |
| mRNA | 19836 | GENCODE V27 |
| snoRNA | 943 | GENCODE V27 |
| snRNA | 1900 | GENCODE V27 |
| srpRNA | 680 | GENCODE V27 |
| tRNA | 649 | GENCODE V27 |
| tucpRNA | 3734 | GENCODE V27 |
| Y\_RNA | 756 | GENCODE V27 |
| circRNA | 140527 | circBase |
| repeats | - | UCSC Genome Browser \(rmsk\) |
| promoter | - | ChromHMM tracks from 9 cell lines from UCSC Genome Browser |
| enhancer | - | ChromHMM tracks from 9 cell lines from UCSC Genome Browser |

spike-in is a special type of genome annotation that should be provided by the user if spike-in sequences are used. 

The paths of the bowtie2 index files:

| Type | FASTA file | bowtie2 index file |
| :--- | :--- | :--- |
| spike-in | `${genome_dir}/fasta/spikein_small.fa` | `${genome_dir}/index/bowtie2/spikein` |
| rRNA | `${genome_dir}/fasta/rRNA.fa` | `${genome_dir}/index/bowtie2/rRNA` |
| miRNA | `${genome_dir}/fasta/miRNA.fa` | `${genome_dir}/rsem_index/bowtie2/miRNA` |
| piRNA | `${genome_dir}/fasta/piRNA.fa` | `${genome_dir}/rsem_index/bowtie2/piRNA` |
| lncRNA | `${genome_dir}/fasta/lncRNA.fa` | `${genome_dir}/rsem_index/bowtie2/lncRNA` |
| mRNA | `${genome_dir}/fasta/mRNA.fa` | `${genome_dir}/rsem_index/bowtie2/mRNA` |
| snoRNA | `${genome_dir}/fasta/snoRNA.fa` | `${genome_dir}/rsem_index/bowtie2/snoRNA` |
| snRNA | `${genome_dir}/fasta/snRNA.fa` | `${genome_dir}/rsem_index/bowtie2/snRNA` |
| srpRNA | `${genome_dir}/fasta/srpRNA.fa` | `${genome_dir}/rsem_index/bowtie2/srpRNA` |
| tRNA | `${genome_dir}/fasta/tRNA.fa` | `${genome_dir}/rsem_index/bowtie2/tRNA` |
| tucpRNA | `${genome_dir}/fasta/tucpRNA.fa` | `${genome_dir}/rsem_index/bowtie2/tucpRNA` |
| Y_RNA | `${genome_dir}/fasta/Y_RNA.fa` | `${genome_dir}/rsem_index/bowtie2/Y_RNA` |
| circRNA | `${genome_dir}/fasta/circRNA.fa` | `${genome_dir}/rsem_index/bowtie2/circRNA` |

**Note**: `${genome_dir}` is the root directory of genome annotation files.

### Build bowtie2 index for spike-in sequences

If your samples contain spike-in sequences, you should first prepare a FASTA file of your spike-in sequences and copy it to `${genome_dir}/fasta/spikein_small.fa`. 
Then create an index file (`${genome_dir}/fasta/spikein_small.fai`) by the following command:

```bash
samtools faidx ${genome_dir}/fasta/spikein_small.fa
```

Run the following commands to build bowtie2 index files for spike-in sequences:

```bash
cut -f1,2 ${genome_dir}/fasta/spikein_small.fa.fai > ${genome_dir}/chrom_sizes/spikein_small
{
    echo -e 'chrom\tstart\tend\tname\tscore\tstrand\tgene_id\ttranscript_id\tgene_name\ttranscript_name\tgene_type\ttranscript_type\tsource'
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,0,$2,$1,0,"+",$1,$1,$1,$1,"spikein","spikein","spikein"}' ${genome_dir}/fasta/spikein_small.fa.fai
} > ${genome_dir}/transcript_table/spikein_small.txt
bowtie2-build ${genome_dir}/fasta/spikein_small.fa ${genome_dir}/index/bowtie2/spikein_small
```
