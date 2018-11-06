#需要更改第一行代码和输入或输出文件的路径
path0=<your path to>/proj_exRNA;
for i in `ls $path0/output/01.trim/cutadapt`;
do j=${i%.cutadapt.fastq};
echo $j;
cleanN=`cat $path0/output/01.trim/cutadapt/$j.cutadapt.fastq | wc -l`
let cleanN=$cleanN/4;
rRNAN=`samtools view $path0/output/02.mapping/1.no_rRNA/sam/$j.rRNA.sam | wc -l`
keptN=`cat $path0/output/02.mapping/1.no_rRNA/fastq/$j.no_rRNA.fq | wc -l`
let keptN=$keptN/4;
miRNAN=`samtools view $path0/output/02.mapping/2.no_miRNA/sam/$j.miRNA.sam | wc -l`
de_miRNAN=`cat $path0/output/02.mapping/2.no_miRNA/fastq/$j.miRNA.unAligned.fastq | wc -l`
let de_miRNAN=$de_miRNAN/4;
piRNAN=`samtools view $path0/output/02.mapping/3.no_piRNA/sam/$j.piRNA.sam | wc -l`
de_piRNAN=`cat $path0/output/02.mapping/3.no_piRNA/fastq/$j.piRNA.unAligned.fastq | wc -l`
let de_piRNAN=$de_piRNAN/4;
Y_RNAN=`samtools view $path0/output/02.mapping/4.no_Y_RNA/sam/$j.Y_RNA.sam | wc -l`
de_Y_RNAN=`cat $path0/output/02.mapping/4.no_Y_RNA/fastq/$j.Y_RNA.unAligned.fastq | wc -l`
let de_Y_RNAN=$de_Y_RNAN/4;
srpRNAN=`samtools view $path0/output/02.mapping/5.no_srpRNA/sam/$j.srpRNA.sam | wc -l`
de_srpRNAN=`cat $path0/output/02.mapping/5.no_srpRNA/fastq/$j.srpRNA.unAligned.fastq | wc -l`
let de_srpRNAN=$de_srpRNAN/4;
tRNAN=`samtools view $path0/output/02.mapping/6.no_tRNA/sam/$j.tRNA.sam | wc -l`
de_tRNAN=`cat $path0/output/02.mapping/6.no_tRNA/fastq/$j.tRNA.unAligned.fastq | wc -l`
let de_tRNAN=$de_tRNAN/4;
snRNAN=`samtools view $path0/output/02.mapping/7.no_snRNA/sam/$j.snRNA.sam | wc -l`
de_snRNAN=`cat $path0/output/02.mapping/7.no_snRNA/fastq/$j.snRNA.unAligned.fastq | wc -l`
let de_snRNAN=$de_snRNAN/4;
snoRNAN=`samtools view $path0/output/02.mapping/8.no_snoRNA/sam/$j.snoRNA.sam | wc -l`
de_snoRNAN=`cat $path0/output/02.mapping/8.no_snoRNA/fastq/$j.snoRNA.unAligned.fastq | wc -l`
let de_snoRNAN=$de_snoRNAN/4;
lncRNAN=`samtools view $path0/output/02.mapping/9.no_lncRNA/sam/$j.lncRNA.sam | wc -l`
de_lncRNAN=`cat $path0/output/02.mapping/9.no_lncRNA/fastq/$j.lncRNA.unAligned.fastq | wc -l`
let de_lncRNAN=$de_lncRNAN/4;
mRNAN=`samtools view $path0/output/02.mapping/10.no_mRNA/sam/$j.mRNA.sam | wc -l`
de_mRNAN=`cat $path0/output/02.mapping/10.no_mRNA/fastq/$j.mRNA.unAligned.fastq | wc -l`
let de_mRNAN=$de_mRNAN/4;
tucpN=`samtools view $path0/output/02.mapping/11.no_tucp/sam/$j.tucp.sam | wc -l`
de_tucpN=`cat $path0/output/02.mapping/11.no_tucp/fastq/$j.tucp.unAligned.fastq | wc -l`
let de_tucpN=$de_tucpN/4;
hg38otherN=`samtools view $path0/output/02.mapping/12.no_hg38other/sam/$j.hg38other.sam | wc -l`
de_hg38otherN=`cat $path0/output/02.mapping/12.hg38other/fastq/$j.hg38other.unAligned.fastq | wc -l`
let de_hg38otherN=$de_hg38otherN/4;
echo -e "$j\t$cleanN\t$rRNAN\t$keptN\t$miRNAN\t$de_miRNAN\t$piRNAN\t$de_piRNAN\t$Y_RNAN\t$de_Y_RNAN\t$srpRNAN\t$de_srpRNAN\t$tRNAN\t$de_tRNAN\t$snRNAN\t$de_snRNAN\t$snoRNAN\t$de_snoRNAN\t$lncRNAN\t$de_lncRNAN\t$mRNAN\t$de_mRNAN\t$tucpN\t$de_tucpN\t$hg38otherN\t$de_hg38otherN\t" | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18,$19,$20,$21,$22,$23,$24,$25,$26}' >> $path0/stat/proj_exRNA.readsN.stat.txt
done;
sed -i -e "1i cleanN\trRNAN\tkeptN\tmiRNAN\tde_miRNAN\tpiRNAN\tde_piRNAN\tY_RNAN\tde_Y_RNAN\tsrpRNAN\tde_srpRNAN\ttRNAN\tde_tRNAN\tsnRNAN\tde_snRNAN\tsnoRNAN\tde_snoRNAN\tlncRNAN\tde_lncRNAN\tmRNAN\tde_mRNAN\ttucpN\tde_tucpN\thg38otherN\tde_hg38otherN" $path0/stat/proj_exRNA.readsN.stat.txt
###比较sam/bam文件差别###

