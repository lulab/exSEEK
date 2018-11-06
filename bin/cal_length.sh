#需要更改第一行代码和输入或输出文件的路径
path0=<your path to>/proj_exRNA;
l=1;
for i in `ls $path0/data/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
for k in rRNA miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp;
do echo $k;
samtools view $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.bam | awk 'BEGIN{FS=OFS="\t"}{print length($10)}' | sort -n | uniq -c | awk 'BEGIN{FS=" "; OFS="\t"}{print $2,$1}' | sort -nk1,1 | sed -e "s/^/$j\t$k\t/g" >> $path0/stat/$j.lengthN.stat.tsv
let "l++";
done;
l=12;
k='hg38other';
echo hg38other;
samtools view $path0/output/02.mapping/$l.no_hg38other/bam/$j.hg38other.bam | awk 'BEGIN{FS=OFS="\t"}{print length($10)}' | sort -n | uniq -c | awk 'BEGIN{FS=" "; OFS="\t"}{print $2,$1}' | sort -nk1,1 | sed -e "s/^/$j\t$k\t/g" >> $path0/stat/$j.lengthN.stat.tsv
l=1;
sed -i -e "1i sample\ttype\tlen\tnum" $path0/stat/$j.lengthN.stat.tsv
done;

