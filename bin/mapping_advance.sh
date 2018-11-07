########cutadpater########
path0=/home/xieyufeng/proj_exRNA;
cd $path0/output/02.mapping;
rm -rf *;
cd $path0/output/01.trim
mkdir short cutadapt
for i in `ls $path0/raw_data`;
do j=${i%.*};
echo $j;
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -m 16 --trim-n --too-short-output=$path0/output/01.trim/short/$j.tooShort.fastq -o $path0/output/01.trim/cutadapt/$j.cutadapt.fastq $path0/raw_data/$j.fastq > {log.log}
done
########rm rRNA########
path0=/home/xieyufeng/proj_exRNA;
mkdir $path0/output/02.mapping/no_rRNA
cd $path0/output/02.mapping/no_rRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/01.trim/cutadapt`;
do j=${i%.cutadapt.fastq};
echo $j;
bowtie2 -p 4 --norc --sensitive-local --no-unal --un $path0/output/02.mapping/no_rRNA/fastq/$j.no_rRNA.fq -x $path0/RNA_index/rRNA $path0/output/01.trim/cutadapt/$j.cutadapt.fastq -S $path0/output/02.mapping/no_rRNA/sam/$j.rRNA.sam
rsem-tbam2gbam $path0/RNA_index/rRNA $path0/output/02.mapping/no_rRNA/sam/$j.rRNA.sam $path0/output/02.mapping/no_rRNA/rsem_bam/$j.rRNA.rsem.clean.bam
done
########rm miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp########
arr=("miRNA" "piRNA" "Y_RNA" "srpRNA" "tRNA" "snRNA" "snoRNA" "lncRNA" "mRNA" "tucp")
for k in {1..10};do
mkdir $path0/output/02.mapping/no_$k
cd $path0/output/02.mapping/no_$k
mkdir fastq sam rsem_bam
for i in `ls $path0/output/01.trim/cutadapt`;
do j=${i%.cutadapt.fastq};
echo $j;
let "l=k-1"
if $l==0;then
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/no_$arr[k]/fastq/$j.$arr[k].unAligned.fastq -x $path0/RNA_index/$arr[k] $path0/output/02.mapping/no_rRNA/fastq/$j.rRNA.unAligned.fastq -S $path0/output/02.mapping/no_$arr[k]/sam/$j.$arr[k].sam
rsem-tbam2gbam $path0/RNA_index/$arr[k] $path0/output/02.mapping/no_$arr[k]/sam/$j.$arr[k].sam $path0/output/02.mapping/no_$arr[k]/rsem_bam/$j.$arr[k].rsem.clean.bam
else
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/no_$arr[k]/fastq/$j.$arr[k].unAligned.fastq -x $path0/RNA_index/$arr[k] $path0/output/02.mapping/no_$arr[l]/fastq/$j.rRNA.unAligned.fastq -S $path0/output/02.mapping/no_$arr[k]/sam/$j.$arr[k].sam
rsem-tbam2gbam $path0/RNA_index/$arr[k] $path0/output/02.mapping/no_$arr[k]/sam/$j.$arr[k].sam $path0/output/02.mapping/no_$arr[k]/rsem_bam/$j.$arr[k].rsem.clean.bam
fi
done

########hg38other########
mkdir $path0/output/02.mapping/hg38other
cd $path0/output/02.mapping/hg38other
mkdir fastq sam bam
for i in `ls $path0/output/01.trim/cutadapt`;
do j=${i%.cutadapt.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/hg38other/fastq/$j.hg38other.unAligned.fastq -x $path0/hg38_index/GRCh38_p10 $path0/output/02.mapping/no_$arr[k]/fastq/$j.$arr[k].unAligned.fastq -S $path0/output/02.mapping/hg38other/sam/$j.hg38other.sam
samtools view -S -b $path0/output/02.mapping/hg38other/sam/$j.hg38other.sam > $path0/output/02.mapping/hg38other/bam/$j.hg38other.bam
done