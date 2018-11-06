path0=/home/xieyufeng/proj_exRNA
## 1 ##
# Build bedGraph files
# use BEDtools and UCSC Kent Utilities
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
l=1;
for k in rRNA miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp;
do echo $k;
samtools sort $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.bam > $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.bam
samtools index $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.bam
bedtools genomecov -ibam $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.bam -bga -split | LC_ALL=C sort -k1,1 -k2,2n > $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.bedGraph
bedGraphToBigWig $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.bedGraph /BioII/lulab_b/shared/genomes/human_hg38/sequence/hg38.chrom.sizes $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.bw
let "l++";
done;
l=12;
samtools sort $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.bam > $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.bam
samtools index $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.bam
bedtools genomecov -ibam $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.bam -bga -split | LC_ALL=C sort -k1,1 -k2,2n > $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.bedGraph
bedGraphToBigWig $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.bedGraph /BioII/lulab_b/shared/genomes/human_hg38/sequence/hg38.chrom.sizes $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.bw
done;

# use homer
mkdir $path0/output/03.tags/
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
l=1;
for k in rRNA miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp;
do echo $k;
mkdir $path0/output/03.tags/$j/$k
makeTagDirectory $path0/output/03.tags/$j/$k $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.bam
# Make bedGraph visualization files for each tag directory
makeUCSCfile $path0/output/03.tags/$j/$k -fragLength given -o auto
let "l++";
done;
l=12;
mkdir $path0/output/03.tags/$j/hg38other
makeTagDirectory $path0/output/03.tags/$j/hg38other $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.bam
makeUCSCfile $path0/output/03.tags/$j/hg38other -fragLength given -o auto
done;

## 2 ##
## calculate raw counts/rpkm/rpm
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
l=1;
for k in rRNA miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp;
do echo $k;
samtools view -h $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.bam > $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.sam
let "l++";
done;
l=12;
samtools view -h $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.bam > $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.sam
done;

# merge alignments of RNA types into one single sam/bam file;
mkdir $path0/output/02.mapping/merged
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
samtools merge $path0/output/02.mapping/merged/$j.merged.bam `find $path0/output/02.mapping/* -name "*rsem.clean.sorted.bam"` $path0/output/02.mapping/12.hg38other/bam/$j.hg38other.sorted.bam
samtools view -h $path0/output/02.mapping/merged/$j.merged.bam > $path0/output/02.mapping/merged/$j.merged.sam

makeTagDirectory $path0/output/03.tags/$j/merged $path0/output/02.mapping/merged/$j.merged.bam
done;



# count for all miRNAs
mkdir $path0/output/04.counts
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
mkdir $path0/output/04.counts/$j
# HTSeq
#htseq-count -m intersection-strict --idattr=ID --type=miRNA_primary_transcript $path0/02.mapping/$j/miRNA/$j.miRNA.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gff > $path0/04.counts/$j/$j.miRNA.htseq.ct
# featureCounts
#featureCounts -t miRNA_primary_transcript -g ID -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gff -o $path0/04.counts/$j/$j.miRNA.featureCounts.counts $path0/02.mapping/$j/miRNA/$j.miRNA.rsem.clean.sorted.bam

l=1;
for RNA_type in rRNA miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp;
do echo $RNA_type;
# featurecounts
featureCounts  -t exon -g transcript_id -s 1 -a $path0/anno/gtf/$RNA_type.* -o $path0/output/04.counts/$j/$j.$RNA_type.featureCounts.counts $path0/output/02.mapping/$l.no_$RNA_type/rsem_bam/$j.$RNA_type.rsem.clean.bam
# raw counts
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/$RNA_type.* hg38 -count exons -d $path0/output/03.tags/$j/$RNA_type/ -gid -noadj > $path0/output/04.counts/$j/$j.$RNA_type.homer.ct 
# calculate rpkm
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/$RNA_type.* hg38 -count exons -d $path0/output/03.tags/$j/merged/ -gid -rpkm > $path0/output/04.counts/$j/$j.$RNA_type.homer.rpkm
# htseq-count
htseq-count -f bam -m intersection-strict $path0/output/02.mapping/$l.no_$RNA_type/rsem_bam/$j.$RNA_type.rsem.clean.bam $path0/anno/gff/$RNA_type.* > $path0/output/04.counts/$j/$j.$RNA_type.htseq.ct
let "l++";
done
# calculate rpm/cpm
#analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.miRNA.homer.rpm
done;




## 3 ##
## merge homer expression matrix 
for k in rRNA miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp;
do echo $k
analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/$k.* hg38 \
-d $path0/output/03.tags/Sample_N1/$k/ $path0/output/03.tags/Sample_N13/$k/ $path0/output/03.tags/Sample_N19/$k/ \
$path0/output/03.tags/Sample_N25/$k/ $path0/output/03.tags/Sample_N7/$k/  \
-gid -noadj > $path0/output/04.counts/proj_exRNA.$k.homer.ct.tsv
cut -f 1,9- $path0/output/04.counts/proj_exRNA.$k.homer.ct.tsv > $path0/output/04.counts/proj_exRNA.$k.homer.ct.mx

analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/anno/gtf/$k.* hg38 \
-d $path0/output/03.tags/Sample_N1/merged/ $path0/output/03.tags/Sample_N13/merged/ $path0/output/03.tags/Sample_N19/merged/ \
$path0/output/03.tags/Sample_N25/merged/ $path0/output/03.tags/Sample_N7/merged/  \
-gid -rpkm > $path0/output/04.counts/proj_exRNA.$k.homer.rpkm.tsv

analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/$k.gtf hg38 \
-d $path0/03.tags/NC_1/merged/ $path0/03.tags/NC_2/merged/ $path0/03.tags/NC_3/merged/ \
$path0/03.tags/BeforeSurgery_1/merged/ $path0/03.tags/BeforeSurgery_2/merged/ $path0/03.tags/BeforeSurgery_3/merged/ \
$path0/03.tags/AfterSurgery_1/merged/ $path0/03.tags/AfterSurgery_2/merged/ $path0/03.tags/AfterSurgery_3/merged/ \
-gid -norm 1e7 > $path0/04.counts/hcc_example.$k.homer.rpm.tsv
done;

## merge htseq expression matrix 
path0=/home/xieyufeng/proj_exRNA

for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;
head -n -5 $path0/output/04.counts/${j}/${j}.${k}.htseq.ct | cut -f 2 | sed -e "1i ${j}" | sed 's/-/_/g' > $path0/output/tmp/${j}.${k}.htseq.ct.tmp
done;
done;

for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;
head -n -5 $path0/output/04.counts/Sample_N13/Sample_N13.${k}.htseq.ct | cut -f 1 | sed -e "s/^/${k}_/g" | sed -e "1i geneID" > $path0/output/tmp/$k.htseq.header
done;
mkdir $path0/output/05.matrix
for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;
paste $path0/output/tmp/$k.htseq.header $path0/output/tmp/*.${k}.htseq.ct.tmp > $path0/output/05.matrix/proj_exRNA.htseq.$k.mx
sed -i "1s/......\t//1" $path0/output/05.matrix/proj_exRNA.htseq.$k.mx
done;

### merge all RNA-type
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;
head -n -5 $path0/output/04.counts/${j}/${j}.${k}.htseq.ct | cut -f 2 > $path0/output/tmp/${j}.${k}.htseq.ct.tmp
done;
cat $path0/output/tmp/${j}.miRNA.htseq.ct.tmp $path0/output/tmp/${j}.piRNA.htseq.ct.tmp $path0/output/tmp/${j}.Y_RNA.htseq.ct.tmp $path0/output/tmp/${j}.snRNA.htseq.ct.tmp $path0/output/tmp/${j}.snoRNA.htseq.ct.tmp $path0/output/tmp/${j}.srpRNA.htseq.ct.tmp $path0/output/tmp/${j}.tRNA.htseq.ct.tmp $path0/output/tmp/${j}.lncRNA.htseq.ct.tmp $path0/output/tmp/${j}.mRNA.htseq.ct.tmp $path0/output/tmp/${j}.tucp.htseq.ct.tmp | sed -e "1i ${j}" | sed 's/-/_/g' > $path0/output/tmp/${j}.merged.htseq.ct.tmp
cat $path0/output/tmp/miRNA.htseq.header $path0/output/tmp/piRNA.htseq.header $path0/output/tmp/Y_RNA.htseq.header $path0/output/tmp/snRNA.htseq.header $path0/output/tmp/snoRNA.htseq.header $path0/output/tmp/srpRNA.htseq.header $path0/output/tmp/tRNA.htseq.header $path0/output/tmp/lncRNA.htseq.header $path0/output/tmp/mRNA.htseq.header $path0/output/tmp/tucp.htseq.header | grep -v "geneID" | sed -e "1i geneID" > $path0/output/tmp/merged.htseq.header
done;

paste $path0/output/tmp/merged.htseq.header $path0/output/tmp/*.merged.htseq.ct.tmp > $path0/output/05.matrix/proj_exRNA.htseq.merged.mx
sed -i "1s/......\t//1" $path0/output/05.matrix/proj_exRNA.htseq.merged.mx



path0=/home/xieyufeng/proj_exRNA
## merge featurecounts expression matrix for each RNA type
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;

head -n -5 $path0/output/04.counts/${j}/${j}.${k}.featureCounts.counts | cut -f 7 | sed -e "1i ${j}" | sed 's/-/_/g' > $path0/output/tmp/${j}.${k}.featureCounts.counts.tmp
sed -i '1,3d' $path0/output/tmp/${j}.${k}.featureCounts.counts.tmp
done;
done;

for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;

head -n -5 $path0/output/04.counts/Sample_N13/Sample_N13.${k}.featureCounts.counts | cut -f 1 > $path0/output/tmp/$k.featureCounts.counts.header
sed -i '1,2d' $path0/output/tmp/$k.featureCounts.counts.header
done;
mkdir $path0/output/05.matrix
for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;
paste $path0/output/tmp/$k.featureCounts.counts.header $path0/output/tmp/*.${k}.featureCounts.counts.tmp > $path0/output/05.matrix/proj_exRNA.featureCounts.counts.$k.mx
#sed -i "1s/......\t//1" $path0/output/05.matrix/proj_exRNA.featureCounts.counts.$k.mx
done;

## merge all RNA-type
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;
head -n -5 $path0/output/04.counts/${j}/${j}.${k}.featureCounts.counts | cut -f 7 | sed -e "1i ${j}" | sed 's/-/_/g'> $path0/output/tmp/${j}.${k}.featureCounts.counts.tmp
sed -i '1,3d' $path0/output/tmp/${j}.${k}.featureCounts.counts.tmp
done;
cat $path0/output/tmp/${j}.miRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.piRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.Y_RNA.featureCounts.counts.tmp $path0/output/tmp/${j}.snRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.snoRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.srpRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.tRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.lncRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.mRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.tucp.featureCounts.counts.tmp | sed -e "1i ${j}" | sed 's/-/_/g' > $path0/output/tmp/${j}.merged.featureCounts.counts.tmp
cat $path0/output/tmp/miRNA.featureCounts.counts.header $path0/output/tmp/piRNA.featureCounts.counts.header $path0/output/tmp/Y_RNA.featureCounts.counts.header $path0/output/tmp/snRNA.featureCounts.counts.header $path0/output/tmp/snoRNA.featureCounts.counts.header $path0/output/tmp/srpRNA.featureCounts.counts.header $path0/output/tmp/tRNA.featureCounts.counts.header $path0/output/tmp/lncRNA.featureCounts.counts.header $path0/output/tmp/mRNA.featureCounts.counts.header $path0/output/tmp/tucp.featureCounts.counts.header | grep -v "geneID" | sed -e "1i geneID" > $path0/output/tmp/merged.featureCounts.counts.header
done;

paste $path0/output/tmp/merged.featureCounts.counts.header $path0/output/tmp/*.merged.featureCounts.counts.tmp > $path0/output/05.matrix/proj_exRNA.featureCounts.counts.merged.mx

#sed -i "1s/......\t//1" $path0/output/05.matrix/proj_exRNA.featureCounts.counts.merged.mx
#sed -i '2,3d' $path0/output/05.matrix/proj_exRNA.featureCounts.counts.merged.mx
#sed -e "s/^/${k}_/g" | sed -e "1i geneID" $path0/output/05.matrix/proj_exRNA.featureCounts.counts.merged.mx

