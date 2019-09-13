# Genome and Annotations

## Annotation summary table

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

## Genome and annotation files

| File | Description |
| :--- | :--- |
| `fasta/genome.fa` | genome sequence |
| `fasta/circRNA.fa` | junction sequence in circBase |
| `fasta/rRNA.fa` | rRNA sequences in NCBI RefSeq |
| `fasta/miRNA.fa` | miRNA hairpin \(precursor\) sequences in miRBase |
| `fasta/piRNA.fa` | piRNA sequences in piRNABank |
| `fasta/${rna_type}.fa` | longest isoform for each gene extracted from GENCODE annotations |
| `gtf_by_biotype/${rna_type}.gtf` | separate GTF files for each RNA type |
| `gtf/gencode.gtf` | GENCODE GTF file |
| `gtf/mitranscriptome.gtf` | Mitranscriptome GTF file |
| `gtf/long_RNA.gtf` | GTF file of Long RNA \(GENCODE + Mitranscriptome - miRNA\) |
| `gtf/piRNABank.gtf` | piRNA GTF file from piRNABank |
| `gtf/gencode_tRNA.gtf` | GTF file of tRNA from GENCODE |
| `transcript_table/all.txt` | table of transcript information \(gene\_id, transcript\_id\) |
| `rsem_index/bowtie2/${rna_type}` | RSEM index files for each RNA type \(built using the longest transcripts\) |
| `rsem_index/bowtie2/${rna_type}.transcripts.fa` | sequence for each RNA type \(longest transcripts\) |
| `gtf_longest_transcript/${rna_type}.gtf` | GTF files for the longest isoforms from GENCODE and Mitranscriptome |
| `bed/*.bed` | transcript in BED12 format extracted from GTF files in \`gtf/\*.gtf |
| `index/bowtie2/${rna_type}` | STAR index for transcripts |
| `index/star/${rna_type}` | STAR index for transcripts |
| `long_index/star/` | STAR index including splicing junctions of long RNA |

## Generate the genome and annotation files

### Create genome directory

```bash
[ -d "genome/hg38/source" ] || mkdir -p "genome/hg38/source"
```

### Chromosome ID conversion table

* Column 1: UCSC chromosome ID
* Column 2: RefSeq chromosome ID

```bash
mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N -e "SELECT * FROM ucscToRefSeq;" hg38 | cut -f1,4 > genome/hg38/source/ucscToRefSeq.txt
```

### Download Gene annotation \(NCBI\)

```bash
# NCBI Human Release 109
wget -P genome/hg38/source ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/GFF/ref_GRCh38.p12_top_level.gff3.gz
[ -d genome/hg38/gff3 ] || mkdir -p genome/hg38/gff3
awk 'BEGIN{OFS="\t";FS="\t"} NR==FNR{c[$2]=$1;next} !/^#/{chrom=c[$1]; if(length(chrom) > 0) print c[$1],$2,$3,$4,$5,$6,$7,$8,$9}' \
    genome/hg38/source/ucscToRefSeq.txt <(zcat genome/hg38/source/ref_GRCh38.p12_top_level.gff3.gz) \
    > genome/hg38/gff3/refseq.gff3
gffread --bed -o genome/hg38/bed/ncbi.bed genome/hg38/gff3/refseq.gff3
wget -O genome/hg38/source/refSeq_rna.fa.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/H_sapiens/RNA/rna.fa.gz
# get rRNA sequence IDs
zgrep 'ribosomal RNA$' genome/hg38/source/refSeq_rna.fa.gz \
    | sed 's/>ref|\(NR_[0-9.]\+\)|.*(\([^)]\+\)).*/\1|\2/' > genome/hg38/source/refSeq_rRNA.ids.txt
# get rRNA sequences
zcat genome/hg38/source/refSeq_rna.fa.gz \
    | awk 'FNR==NR{split($0,a,"|");ids[a[1]]=1;next} 
{if($0 ~ /^>/){split($0,a,"|");if(ids[a[2]] == 1){keep=1; print ">" a[2];}else{keep=0}} else{if(keep == 1){print}}}' \
genome/hg38/source/refSeq_rRNA.ids.txt - > genome/hg38/fasta/rRNA.fa
samtools faidx genome/hg38/fasta/rRNA.fa
# generate transcript table
{
echo -e 'chrom\tstart\tend\tname\tscore\tstrand\tgene_id\ttranscript_id\tgene_name\ttranscript_name\tgene_type\ttranscript_type\tsource'
awk 'BEGIN{OFS="\t";FS="\t"}FNR==NR{split($0,a,"|");gene_name[a[1]]=a[2];next}{print $1,0,$2,a[1],0,"+",$1,$1,gene_name[$1],gene_name[$1],"rRNA","rRNA","RefSeq"}' \
    genome/hg38/source/refSeq_rRNA.ids.txt genome/hg38/fasta/rRNA.fa.fai
} > genome/hg38/transcript_table/rRNA.txt
# get transcript sizes
cut -f1,2 genome/hg38/fasta/rRNA.fa.fai > genome/hg38/chrom_sizes/rRNA
# build STAR index (small genome)
STAR --runMode genomeGenerate --genomeSAindexNbases 5 --genomeDir genome/hg38/index/star/rRNA/ --genomeFastaFiles genome/hg38/fasta/rRNA.fa
```

### Download chain files for CrossMap

```bash
wget -O genome/hg38/source/hg18ToHg38.over.chain.gz http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg38.over.chain.gz
wget -O genome/hg38/source/NCBI36_to_GRCh38.chain.gz https://sourceforge.net/projects/crossmap/files/Ensembl_chain_files/homo_sapiens%28human%29/NCBI36_to_GRCh38.chain.gz
```

### Genome assembly \(UCSC hg38\)

```bash
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gzip -d -c genome/hg38/source/hg38.fa.gz > genome/hg38/fasta/genome.fa
samtools faidx genome/hg38/fasta/genome.fa
```

### ENCODE annotations

```bash
wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gtf.gz
#wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.annotation.gff3.gz
wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.long_noncoding_RNAs.gtf.gz
#wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.long_noncoding_RNAs.gff3.gz
wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.tRNAs.gtf.gz
#wget -P genome/hg38/source ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_27/gencode.v27.tRNAs.gff3.gz
zcat genome/hg38/source/gencode.v27.annotation.gtf.gz > genome/hg38/gtf/gencode.gtf
zcat genome/hg38/source/gencode.v27.long_noncoding_RNAs.gtf.gz > genome/hg38/gtf/gencode_lncRNA.gtf
zcat genome/hg38/source/gencode.v27.tRNAs.gtf.gz \
    | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$2,"exon",$4,$5,$6,$7,$8,$9}' > genome/hg38/gtf/gencode_tRNA.gtf
# Chain file for converting hg19 to hg38
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
```

### Mitranscriptome

```bash
wget -P genome/hg38/source http://mitranscriptome.org/download/mitranscriptome.gtf.tar.gz
tar -C genome/hg38/source --strip-components=1 -zxf genome/hg38/source/mitranscriptome.gtf.tar.gz mitranscriptome.gtf/mitranscriptome.v2.gtf.gz
# convert from hg19 to hg38
zcat genome/hg38/source/mitranscriptome.v2.gtf.gz \
    | CrossMap.py gff genome/hg38/source/hg19ToHg38.over.chain.gz /dev/stdin genome/hg38/source/mitranscriptome.v2.hg38.gtf
# remove invalid transcripts
bin/preprocess.py fix_gtf -i genome/hg38/source/mitranscriptome.v2.hg38.gtf -o genome/hg38/gtf/mitranscriptome.gtf
```

Extract lncRNA and TUCP RNA to separate GTF files:

```bash
grep 'tcat "lncrna"' genome/hg38/gtf/mitranscriptome.gtf > genome/hg38/gtf/mitranscriptome_lncRNA.gtf
# add exon feature
grep 'tcat "tucp"' genome/hg38/gtf/mitranscriptome.gtf \
    | awk 'BEGIN{OFS="\t";FS="\t"}{print;print $1,$2,"exon",$4,$5,$6,$7,$8,$9}' > genome/hg38/gtf/mitranscriptome_tucp.gtf
cp genome/hg38/gtf/mitranscriptome_tucp.gtf genome/hg38/gtf_by_biotype/tucpRNA.gtf
```

### NONCODE

```bash
wget -P genome/hg38/source http://www.noncode.org/datadownload/NONCODEv5_human_hg38_lncRNA.gtf.gz
zcat genome/hg38/source/NONCODEv5_human_hg38_lncRNA.gtf.gz \
    | awk 'BEGIN{FS="\t";OFS="\t"}$7 != "." {print $1,"NONCODE",$3,$4,$5,$6,$7,$8,$9}' > genome/hg38/gtf/noncode.gtf
```

### lncRNAs identified in HCC \(Nature communications 2017\)

```bash
wget -P genome/hg38/source https://media.nature.com/original/nature-assets/ncomms/2017/170213/ncomms14421/extref/ncomms14421-s3.txt
awk 'BEGIN{FS="\t";OFS="\t"}{print $1,"ncomms2017",$3,$4,$5,$6,$7,$8,$9}' genome/hg38/source/ncomms14421-s3.txt > genome/hg38/source/ncomms2017.gtf
CrossMap.py gff genome/hg38/source/hg19ToHg38.over.chain.gz genome/hg38/source/ncomms2017.gtf genome/hg38/source/ncomms2017.hg38.gtf
ln genome/hg38/source/ncomms2017.hg38.gtf genome/hg38/gtf/ncomms2017.gtf
```

### Merge lncRNA \(GENCODE and Mitranscriptome\)

```bash
cat genome/hg38/gtf/gencode_lncRNA.gtf \
    genome/hg38/gtf/mitranscriptome_lncRNA.gtf \
    > genome/hg38/gtf/merged_lncRNA.gtf
cp genome/hg38/gtf/merged_lncRNA.gtf genome/hg38/gtf_by_biotype/lncRNA.gtf
```

### piRBase \(v1.0\)

```bash
wget -O genome/hg38/source/piRBase-hsa-v1.0.bed.gz http://www.regulatoryrna.org/database/piRNA/download/archive/v1.0/bed/piR_hg19_v1.0.bed.gz
zcat genome/hg38/source/piRBase-hsa-v1.0.bed.gz \
    | CrossMap.py bed genome/hg38/source/hg19ToHg38.over.chain.gz /dev/stdin genome/hg38/source/piRBase-hsa-v1.0.hg38.bed
bedToGenePred genome/hg38/source/piRBase-hsa-v1.0.hg38.bed genome/hg38/source/piRBase-hsa-v1.0.hg38.genePred
genePredToGtf -source=piRBase file genome/hg38/source/piRBase-hsa-v1.0.hg38.genePred genome/hg38/source/piRBase-hsa-v1.0.hg38.gtf
ln genome/hg38/source/piRBase-hsa-v1.0.hg38.gtf genome/hg38/gtf/piRBase.gtf
```

### piRBase \(v2.0\)

```bash
wget -O genome/hg38/source/piRBase-hsa-v2.0.bed.gz http://www.regulatoryrna.org/database/piRNA/download/archive/v2.0/bed/hsa.bed.gz
zcat genome/hg38/source/piRBase-hsa-v2.0.bed.gz | bedtools sort > source/piRBase-hsa-v2.0.bed
bedToGenePred source/piRBase-hsa-v2.0.bed source/piRBase-hsa-v2.0.genePred
genePredToGtf -source=piRBase file source/piRBase-hsa-v2.0.genePred source/piRBase-hsa-v2.0.gtf
```

### Long RNA \(GENCODE + Mitranscriptome - miRNA\)

```bash
# Merge GTF files
cat genome/hg38/gtf/gencode.gtf \
    genome/hg38/gtf/mitranscriptome_lncRNA.gtf \
    genome/hg38/gtf/mitranscriptome_tucp.gtf \
    | grep -v 'gene_type "miRNA' \
    > genome/hg38/gtf/long_RNA.gtf
# Get gene lengths
tools/GTFtools_0.6.5/gtftools.py -c 1-22,X,Y,M -l genome/hg38/gene_length/long_RNA genome/hg38/gtf/long_RNA.gtf
# GTF to BED12 format
gffread --bed -o genome/hg38/bed/long_RNA.bed genome/hg38/gtf/long_RNA.gtf
```

**gene\_length/long\_RNA**

* Tab-deliminated text file
* First row: header
* Column 1 \(gene\): gene\_id
* Column 2 \(mean\): mean length of isoforms
* Column 3 \(median\): median length of isoforms
* Column 4 \(longest\_isoform\): length of the longest isoform
* Column 5 \(merged\): merged length of isoforms

### piRNABank \(NCBI36\)

```bash
wget -O genome/hg38/source/ http://pirnabank.ibab.ac.in/downloads/all/human_all.zip
unzip genome/hg38/source/human_all.zip -d genome/hg38/source/
mv genome/hg38/source/human_pir.txt genome/hg38/source/piRNABank.human.txt
# Extract genomic coordinates from piRNABank
awk 'BEGIN{OFS="\t"}
    /^>/{na=split(substr($0,2),a,"|");split(a[na],b,":"); 
    if(b[5]=="Plus"){s="+"} else{s="-"}
    if(a[1]!=name){print "chr" b[2],b[3]-1,b[4],a[1],0,s}
    name=a[1]}' genome/hg38/source/piRNABank.human.txt \
    | bedtools sort > genome/hg38/source/piRNABank.human.bed
awk 'BEGIN{OFS="\t"}
    {if($0 ~ /^>/) {split(substr($0,2),a,"|"); 
        if((a[1] != name)&&(length(seq) > 0)){print ">" name;gsub(/U/,"T",seq);print seq} name=a[1]}
    else{seq=$0}}' genome/hg38/source/piRNABank.human.txt > genome/hg38/source/piRNABank.human.fa
bedToGenePred genome/hg38/source/piRNABank.human.bed genome/hg38/source/piRNABank.human.genePred
genePredToGtf -source=piRNABank file genome/hg38/source/piRNABank.human.genePred stdout \
    | awk '$3=="exon"' > genome/hg38/source/piRNABank.human.gtf

CrossMap.py gff genome/hg38/source/hg18ToHg38.over.chain.gz genome/hg38/source/piRNABank.human.gtf \
    genome/hg38/source/piRNABank.human.hg38.gtf
cp genome/hg38/source/piRNABank.human.hg38.gtf genome/hg38/gtf/piRNABank.gtf
cp genome/hg38/gtf/piRNABank.gtf genome/hg38/gtf_by_biotype/piRNA.gtf
gffread --bed -o genome/hg38/source/piRNABank.human.hg38.bed genome/hg38/source/piRNABank.human.hg38.gtf
bedtools getfasta -s -name -fi genome/hg38/fasta/genome.fa -bed genome/hg38/source/piRNABank.human.hg38.bed -split \
    > genome/hg38/source/piRNABank.human.hg38.fa
```

### miRBase \(Version 22\)

```bash
wget -O genome/hg38/source/miRBase.hsa.gff3 ftp://mirbase.org/pub/mirbase/CURRENT/genomes/hsa.gff3
wget -O genome/hg38/source/miRBase.hairpin.fa.gz ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz
wget -O genome/hg38/source/miRBase.mature.fa.gz ftp://mirbase.org/pub/mirbase/CURRENT/mature.fa.gz

cp genome/hg38/source/miRBase.hsa.gff3 genome/hg38/gtf/miRBase.gff3
# extract human pre-miRNA
zcat genome/hg38/source/miRBase.hairpin.fa.gz \
    | awk '/^>/{if($0 ~ />hsa-/) {keep=1; print $1} else{keep=0}; next}{if(keep==1){gsub(/U/, "T");print}}' \
    > genome/hg38/fasta/miRNA.fa
samtools faidx genome/hg38/fasta/miRNA.fa
# extract human mature miRNA
zcat genome/hg38/source/miRBase.mature.fa.gz \
    | awk '/^>/{if($0 ~ />hsa-/) {keep=1; print $1} else{keep=0}; next}{if(keep==1){gsub(/U/, "T");print}}' \
    > genome/hg38/fasta/mature_miRNA.fa
samtools faidx genome/hg38/fasta/mature_miRNA.fa
# generate transcript table (mature miRNA)
{
echo -e 'chrom\tstart\tend\tname\tscore\tstrand\tgene_id\ttranscript_id\tgene_name\ttranscript_name\tgene_type\ttranscript_type\tsource'
awk 'BEGIN{OFS="\t";FS="\t"}{print $1,0,$2,$1,0,"+",$1,$1,$1,$1,"miRNA","miRNA","miRBase"}' \
    genome/hg38/fasta/miRNA.fa.fai genome/hg38/fasta/mature_miRNA.fa.fai 
} > genome/hg38/transcript_table/miRNA.txt
# get transcript sizes
cut -f1,2 genome/hg38/fasta/miRNA.fa.fai genome/hg38/fasta/mature_miRNA.fa.fai > genome/hg38/chrom_sizes/miRNA
# gff3 to genePred
awk 'BEGIN{OFS="\t";FS="\t";d["miRNA"]="transcript";d["miRNA_primary_transcript"]="primary_transcript"}/^#/{print}!/^#/{$3=d[$3];print $1,$2,$3,$4,$5,$6,$7,$8,$9}' \
    genome/hg38/gff3/miRBase.gff3 > genome/hg38/source/miRBase.fixed.gff3
gff3ToGenePred -useName genome/hg38/source/miRBase.fixed.gff3 genome/hg38/genePred/miRBase.genePred
```

### Spike-in

```bash
cut -f1,2 genome/hg38/fasta/spikein_small.fa.fai > genome/hg38/chrom_sizes/spikein_small
{
    echo -e 'chrom\tstart\tend\tname\tscore\tstrand\tgene_id\ttranscript_id\tgene_name\ttranscript_name\tgene_type\ttranscript_type\tsource'
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,0,$2,$1,0,"+",$1,$1,$1,$1,"spikein","spikein","spikein"}' genome/hg38/fasta/spikein_small.fa.fai
} > genome/hg38/transcript_table/spikein_small.txt
bowtie2-build genome/hg38/fasta/spikein_small.fa genome/hg38/index/bowtie2/spikein_small
STAR --runMode genomeGenerate --genomeSAindexNbases 10 --genomeChrBinNbits 7 \
    --genomeDir genome/hg38/index/star/spikein_long/ --genomeFastaFiles genome/hg38/fasta/spikein_long.fa
```

### UniVec

```bash
wget -O genome/hg38/fasta/univec.fa 'ftp://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec'
samtools faidx genome/hg38/fasta/univec.fa
cut -f1,2 genome/hg38/fasta/univec.fa.fai > genome/hg38/chrom_sizes/univec
{
    echo -e 'chrom\tstart\tend\tname\tscore\tstrand\tgene_id\ttranscript_id\tgene_name\ttranscript_name\tgene_type\ttranscript_type\tsource'
    awk 'BEGIN{OFS="\t";FS="\t"}{print $1,0,$2,$1,0,"+",$1,$1,$1,$1,"univec","univec","univec"}' genome/hg38/fasta/univec.fa.fai
} > genome/hg38/transcript_table/univec.txt
bowtie2-build genome/hg38/fasta/univec.fa genome/hg38/index/bowtie2/univec
STAR --runMode genomeGenerate --genomeSAindexNbases 10 --genomeChrBinNbits 7 --genomeDir genome/hg38/index/star/univec/ --genomeFastaFiles genome/hg38/fasta/univec.fa
```

### Intron

```bash
bin/preprocess.py extract_gene -i genome/hg38/gtf/long_RNA.gtf | bedtools sort > genome/hg38/bed/long_RNA.gene.bed
awk 'BEGIN{OFS="\t";FS="\t"} !/^#/{match($9,/gene_id "([^"]+)"/,a);print $1,$4-1,$5,a[1],0,$7}' genome/hg38/gtf/long_RNA.gtf \
    | bedtools sort > genome/hg38/bed/long_RNA.exon.bed
bedtools subtract -sorted -s -a genome/hg38/bed/long_RNA.gene.bed -b genome/hg38/bed/long_RNA.exon.bed \
    | bedtools sort > genome/hg38/bed/long_RNA.intron.bed
```

### Promoter/enhancer from ChromHMM \(hg19\)

```bash
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmGm12878HMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmH1hescHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHepg2HMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHmecHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHsmmHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmHuvecHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmK562HMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhekHMM.bed.gz
wget -P genome/hg38/source http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeBroadHmm/wgEncodeBroadHmmNhlfHMM.bed.gz

# hg19 => hg38
tracks="wgEncodeBroadHmmGm12878HMM wgEncodeBroadHmmH1hescHMM wgEncodeBroadHmmHepg2HMM
    wgEncodeBroadHmmHmecHMM wgEncodeBroadHmmHsmmHMM wgEncodeBroadHmmHuvecHMM
    wgEncodeBroadHmmK562HMM wgEncodeBroadHmmNhekHMM wgEncodeBroadHmmNhlfHMM"
for track in $tracks;do
    CrossMap.py bed genome/hg38/source/hg18ToHg38.over.chain.gz <(zcat genome/hg38/source/${track}.bed.gz) genome/hg38/source/${track}.hg38.bed
    awk 'BEGIN{OFS="\t";FS="\t"}($4=="1_Active_Promoter")||($4=="2_Weak_Promoter")||($4=="3_Poised_Promoter"){print $1,$2,$3,$4,$5,$6}' \
        genome/hg38/source/${track}.hg38.bed | bedtools sort > genome/hg38/bed/promoter.${track}.bed
    awk 'BEGIN{OFS="\t";FS="\t"}($4=="4_Strong_Enhancer")||($4=="5_Strong_Enhancer")||($4=="6_Weak_Enhancer")||($4=="7_Weak_Enhancer"){print $1,$2,$3,$4,$5,$6}' \
        genome/hg38/source/${track}.hg38.bed | bedtools sort > genome/hg38/bed/enhancer.${track}.bed
done
# merge promoters and enhancers from 9 cell lines
cat $(for track in $tracks;do echo genome/hg38/bed/promoter.${track}.bed;done) \
    | bedtools sort | bedtools merge -d 1 \
    | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,"promoter",0,"."}' > genome/hg38/bed/promoter.merged.bed
cat $(for track in $tracks;do echo genome/hg38/bed/enhancer.${track}.bed;done) \
    | bedtools sort | bedtools merge -d 1 \
    | awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,"enhancer",0,"."}' > genome/hg38/bed/enhancer.merged.bed
# convert unstranded bed to stranded bed (duplicate each record into two records with "+" and "-")
awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,$4,$5,"+";print $1,$2,$3,$4,$5,"-"}' \
    genome/hg38/bed/promoter.merged.bed > genome/hg38/bed/promoter.stranded.bed
awk 'BEGIN{OFS="\t";FS="\t"}{print $1,$2,$3,$4,$5,"+";print $1,$2,$3,$4,$5,"-"}' \
    genome/hg38/bed/enhancer.merged.bed > genome/hg38/bed/enhancer.stranded.bed
```

### Repeats

UCSC GenomeBrowser -&gt; Tools -&gt; Table Browser

* assembly: GRCh38/hg38
* group: repeats
* track: RepeatMasker
* table: rmsk

Dowload to: genome/hg38/source/rmsk.bed.gz

```bash
gunzip -c genome/hg38/source/rmsk.bed.gz | bedtools sort > genome/hg38/bed/rmsk.bed
```

### circRNA database \(circBase\)

```bash
wget -O genome/hg38/source/circbase.hg19.fa.gz http://www.circbase.org/download/human_hg19_circRNAs_putative_spliced_sequence.fa.gz
zcat genome/hg38/source/circbase.hg19.fa.gz | bin/preprocess.py extract_circrna_junction -s 150 -o genome/hg38/fasta/circRNA.fa
samtools faidx genome/hg38/fasta/circRNA.fa
STAR --runMode genomeGenerate --genomeSAindexNbases 10 --genomeChrBinNbits 7 --genomeDir genome/hg38/index/star/circRNA/ --genomeFastaFiles genome/hg38/fasta/circRNA.fa
```

### Create pseudo-genome for IGV

```bash
bin/preprocess.py create_pseudo_genome \
    -i genome/hg38/fasta/circRNA.fa \
    --circular-rna \
    --output-fasta genome/hg38/fasta/pseudo_genome.circRNA.fa \
    --output-annotation genome/hg38/bed/pseudo_genome.circRNA.bed \
    --output-chrom-sizes genome/hg38/chrom_sizes/pseudo_genome.circRNA \
    --output-cytoband genome/hg38/cytoband/pseudo_genome.circRNA.txt
bedtools sort -i genome/hg38/bed/pseudo_genome.circRNA.bed
```

### Merge transcript table

```bash
{
echo -e 'chrom\tstart\tend\tname\tscore\tstrand\tgene_id\ttranscript_id\tgene_name\ttranscript_name\tgene_type\ttranscript_type\tsource'
for rna_type in rRNA lncRNA miRNA mRNA piRNA snoRNA snRNA srpRNA tRNA tucpRNA Y_RNA;do
    sed '1 d' genome/hg38/transcript_table/${rna_type}.txt
done
} > genome/hg38/transcript_table/all.txt
```

