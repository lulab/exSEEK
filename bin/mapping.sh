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
mkdir $path0/output/02.mapping/1.no_rRNA
cd $path0/output/02.mapping/1.no_rRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/01.trim/cutadapt`;
do j=${i%.cutadapt.fastq};
echo $j;
bowtie2 -p 4 --norc --sensitive-local --no-unal --un $path0/output/02.mapping/1.no_rRNA/fastq/$j.no_rRNA.fq -x $path0/RNA_index/rRNA $path0/output/01.trim/cutadapt/$j.cutadapt.fastq -S $path0/output/02.mapping/1.no_rRNA/sam/$j.rRNA.sam
rsem-tbam2gbam $path0/RNA_index/rRNA $path0/output/02.mapping/1.no_rRNA/sam/$j.rRNA.sam $path0/output/02.mapping/1.no_rRNA/rsem_bam/$j.rRNA.rsem.clean.bam
done
########rm miRNA########
mkdir $path0/output/02.mapping/2.no_miRNA
cd $path0/output/02.mapping/2.no_miRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/1.no_rRNA/fastq/`;
do j=${i%.no_rRNA.fq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/2.no_miRNA/fastq/$j.miRNA.unAligned.fastq -x $path0/RNA_index/miRNA $path0/output/02.mapping/1.no_rRNA/fastq/$j.no_rRNA.fq -S $path0/output/02.mapping/2.no_miRNA/sam/$j.miRNA.sam
rsem-tbam2gbam $path0/RNA_index/miRNA $path0/output/02.mapping/2.no_miRNA/sam/$j.miRNA.sam $path0/output/02.mapping/2.no_miRNA/rsem_bam/$j.miRNA.rsem.clean.bam
done
########rm piRNA########
mkdir $path0/output/02.mapping/3.no_piRNA
cd $path0/output/02.mapping/3.no_piRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/2.no_miRNA/fastq/`;
do j=${i%.miRNA.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/3.no_piRNA/fastq/$j.piRNA.unAligned.fastq -x $path0/RNA_index/piRNA $path0/output/02.mapping/2.no_miRNA/fastq/$j.miRNA.unAligned.fastq -S $path0/output/02.mapping/3.no_piRNA/sam/$j.piRNA.sam
rsem-tbam2gbam $path0/RNA_index/piRNA $path0/output/02.mapping/3.no_piRNA/sam/$j.piRNA.sam $path0/output/02.mapping/3.no_piRNA/rsem_bam/$j.piRNA.rsem.clean.bam
done
########rm Y_RNA########
mkdir $path0/output/02.mapping/4.no_Y_RNA
cd $path0/output/02.mapping/4.no_Y_RNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/3.no_piRNA/fastq/`;
do j=${i%.piRNA.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/4.no_Y_RNA/fastq/$j.Y_RNA.unAligned.fastq -x $path0/RNA_index/Y_RNA $path0/output/02.mapping/3.no_piRNA/fastq/$j.piRNA.unAligned.fastq -S $path0/output/02.mapping/4.no_Y_RNA/sam/$j.Y_RNA.sam
rsem-tbam2gbam $path0/RNA_index/Y_RNA $path0/output/02.mapping/4.no_Y_RNA/sam/$j.Y_RNA.sam $path0/output/02.mapping/4.no_Y_RNA/rsem_bam/$j.Y_RNA.rsem.clean.bam
done
########rm srpRNA########
mkdir $path0/output/02.mapping/5.no_srpRNA
cd $path0/output/02.mapping/5.no_srpRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/4.no_Y_RNA/fastq/`;
do j=${i%.Y_RNA.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/5.no_srpRNA/fastq/$j.srpRNA.unAligned.fastq -x $path0/RNA_index/srpRNA $path0/output/02.mapping/4.no_Y_RNA/fastq/$j.Y_RNA.unAligned.fastq -S $path0/output/02.mapping/5.no_srpRNA/sam/$j.srpRNA.sam
rsem-tbam2gbam $path0/RNA_index/srpRNA $path0/output/02.mapping/5.no_srpRNA/sam/$j.srpRNA.sam $path0/output/02.mapping/5.no_srpRNA/rsem_bam/$j.srpRNA.rsem.clean.bam
done
########rm tRNA########
mkdir $path0/output/02.mapping/6.no_tRNA
cd $path0/output/02.mapping/6.no_tRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/5.no_srpRNA/fastq/`;
do j=${i%.srpRNA.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/6.no_tRNA/fastq/$j.tRNA.unAligned.fastq -x $path0/RNA_index/tRNA $path0/output/02.mapping/5.no_srpRNA/fastq/$j.srpRNA.unAligned.fastq -S $path0/output/02.mapping/6.no_tRNA/sam/$j.tRNA.sam
rsem-tbam2gbam $path0/RNA_index/tRNA $path0/output/02.mapping/6.no_tRNA/sam/$j.tRNA.sam $path0/output/02.mapping/6.no_tRNA/rsem_bam/$j.tRNA.rsem.clean.bam
done
########rm snRNA########
mkdir $path0/output/02.mapping/7.no_snRNA
cd $path0/output/02.mapping/7.no_snRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/6.no_tRNA/fastq/`;
do j=${i%.tRNA.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/7.no_snRNA/fastq/$j.snRNA.unAligned.fastq -x $path0/RNA_index/snRNA $path0/output/02.mapping/6.no_tRNA/fastq/$j.tRNA.unAligned.fastq -S $path0/output/02.mapping/7.no_snRNA/sam/$j.snRNA.sam
rsem-tbam2gbam $path0/RNA_index/snRNA $path0/output/02.mapping/7.no_snRNA/sam/$j.snRNA.sam $path0/output/02.mapping/7.no_snRNA/rsem_bam/$j.snRNA.rsem.clean.bam
done
########rm snoRNA########
mkdir $path0/output/02.mapping/8.no_snoRNA
cd $path0/output/02.mapping/8.no_snoRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/7.no_snRNA/fastq/`;
do j=${i%.snRNA.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/8.no_snoRNA/fastq/$j.snoRNA.unAligned.fastq -x $path0/RNA_index/snoRNA $path0/output/02.mapping/7.no_snRNA/fastq/$j.snRNA.unAligned.fastq -S $path0/output/02.mapping/8.no_snoRNA/sam/$j.snoRNA.sam
rsem-tbam2gbam $path0/RNA_index/snoRNA $path0/output/02.mapping/8.no_snoRNA/sam/$j.snoRNA.sam $path0/output/02.mapping/8.no_snoRNA/rsem_bam/$j.snoRNA.rsem.clean.bam
done
########rm lncRNA########
mkdir $path0/output/02.mapping/9.no_lncRNA
cd $path0/output/02.mapping/9.no_lncRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/8.no_snoRNA/fastq/`;
do j=${i%.snoRNA.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/9.no_lncRNA/fastq/$j.lncRNA.unAligned.fastq -x $path0/RNA_index/lncRNA $path0/output/02.mapping/8.no_snoRNA/fastq/$j.snoRNA.unAligned.fastq -S $path0/output/02.mapping/9.no_lncRNA/sam/$j.lncRNA.sam
rsem-tbam2gbam $path0/RNA_index/lncRNA $path0/output/02.mapping/9.no_lncRNA/sam/$j.lncRNA.sam $path0/output/02.mapping/9.no_lncRNA/rsem_bam/$j.lncRNA.rsem.clean.bam
done
########rm mRNA########
mkdir $path0/output/02.mapping/10.no_mRNA
cd $path0/output/02.mapping/10.no_mRNA
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/9.no_lncRNA/fastq/`;
do j=${i%.lncRNA.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/10.no_mRNA/fastq/$j.mRNA.unAligned.fastq -x $path0/RNA_index/mRNA $path0/output/02.mapping/9.no_lncRNA/fastq/$j.lncRNA.unAligned.fastq -S $path0/output/02.mapping/10.no_mRNA/sam/$j.mRNA.sam
rsem-tbam2gbam $path0/RNA_index/mRNA $path0/output/02.mapping/10.no_mRNA/sam/$j.mRNA.sam $path0/output/02.mapping/10.no_mRNA/rsem_bam/$j.mRNA.rsem.clean.bam
done
########rm tucpRNA########
mkdir $path0/output/02.mapping/11.no_tucp
cd $path0/output/02.mapping/11.no_tucp
mkdir fastq sam rsem_bam
for i in `ls $path0/output/02.mapping/10.no_mRNA/fastq/`;
do j=${i%.mRNA.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/11.no_tucp/fastq/$j.tucp.unAligned.fastq -x $path0/RNA_index/tucp $path0/output/02.mapping/10.no_mRNA/fastq/$j.mRNA.unAligned.fastq -S $path0/output/02.mapping/11.no_tucp/sam/$j.tucp.sam
rsem-tbam2gbam $path0/RNA_index/tucp $path0/output/02.mapping/11.no_tucp/sam/$j.tucp.sam $path0/output/02.mapping/11.no_tucp/rsem_bam/$j.tucp.rsem.clean.bam
done
########hg38other########
mkdir $path0/output/02.mapping/12.hg38other
cd $path0/output/02.mapping/12.hg38other
mkdir fastq sam bam
for i in `ls $path0/output/02.mapping/11.no_tucp/fastq/`;
do j=${i%.tucp.unAligned.fastq};
echo $j;
bowtie2 -p 4 --sensitive-local --norc --no-unal --un $path0/output/02.mapping/12.hg38other/fastq/$j.hg38other.unAligned.fastq -x $path0/hg38_index/GRCh38_p10 $path0/output/02.mapping/11.no_tucp/fastq/$j.tucp.unAligned.fastq -S $path0/output/02.mapping/12.hg38other/sam/$j.hg38other.sam
samtools view -S -b $path0/output/02.mapping/12.hg38other/sam/$j.hg38other.sam > $path0/output/02.mapping/12.hg38other/bam/$j.hg38other.bam
done