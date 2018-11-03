## 2\) Reads Processing and Mapping 指南   

完成五个样本`Sample_N1, Sample_N7, Sample_N13, Sample_N19, Sample_N25`的mapping和RNA ratio与length的统计工作，其中产生的bam文件供下一步构建expression matrix使用。

**总体流程图**

![](../assets/mapping_pipe%20%281%29.png)

### 2.1\) Data Structure

```text
~/proj_exRNA/
|-- data
    |-- RNA_index      #4.2.2.4 d. 中比对到各种RNA类型的index
    |-- hg38_index     #4.2.2.4 d. 中最后一步所需要的index
    |-- raw_data
|-- stat               #存放最终步骤的信息统计文件
|-- output             #可以根据自己的习惯对output进行整理，以下是按照流程设置output的路径
eg:
    |-- 01.trim        #对应4.2.2.2
        |-- QC1        #对应4.2.2.2 b. step one
        |-- trim       #对应4.2.2.2 b. step two
        |-- QC2        #对应4.2.2.2 b. step three
    |-- 02.mapping     #对应4.2.2.3 c.和4.2.2.4 d.
      |-- 1.no_rRNA
          |-- fastq    #存*.no_rRNA.fq，详见4.2.2.3 c.
          |-- sam      #存*.<rRNA>.sam，详见4.2.2.3 c.
          |-- rsem_bam #将.sam转化为.bam文件，详见4.2.2.3 c.
      |-- 2.no_miRNA   
      |-- ...
      |-- 12.no_hg38other
          |-- fastq    
          |-- sam      
          |-- bam      #.sam转.bam工具不同，文件夹由rsem_bam改名至bam
    |-- 03.tags        #homer构建表达矩阵所需路径，本教程不需要建立此路径
        |-- Sample_N1
            |-- miRNA
            |-- ...
        |-- ...
    |-- 04.counts      #构建表达矩阵
    |-- 05.matrix      #构建表达矩阵
    |-- tmp            #存放中间文件
```

**Inputs**

| **File format** | **Information contained in file** | **File description** |
| :--- | :--- | :--- |
| fastq | **reads** | five samples, GEO link: GSE71008 |

**Outputs**

| **File format** | **Information contained in file** |
| :--- | :--- |
| sam/bam | mapped reads to different kinds of indexes |
| tsv format | stats of RNA ratio and length |

### 2.2\) Running Steps

#### **2.2.1\) 获取数据**

从`/BioII/chenxupeng/student/`上获取基因组数据`hg38`，基因组注释数据`/gtf`，索引文件`/RNA_index`以及原始数据`(fastq files)`到自己的账号下

| data | path |
| :--- | :--- |
| `hg38` | `/BioII/chenxupeng/student/data/hg38_index/GRCh38.p10.genome.fa` |
| `gtf` | `/BioII/chenxupeng/student/data/gtf` |
| `RNA index` | `/BioII/chenxupeng/student/data/RNA_index/` |
| `raw data` | `/BioII/chenxupeng/student/data/raw_data/*.fastq` |

推荐使用`ln`或`cp`命令

#### **2.2.2\) QC-Trim-QC**

这步操作目的主要有两个，一个是检查数据的质量，另一个是减掉接头序列

* **Step one - QC of raw data**

**Input:**

| data type | path |
| :--- | :--- |
| `raw data` | `/BioII/chenxupeng/student/data/raw_data/*.fastq` |

**Software/Parameters:**

`fastqc`

| `options` | function |
| :--- | :--- |
| `-q --quiet` | Supress all progress messages on stdout and only report errors. |
| `-o --outdir` | Create all output files in the specified output directory. |
| `-h --help` | detailed introduction of options |

**Output:**

QC files

* **step two - cut adaptor & trim long read**

**Input:**

| data type | **path** |
| :--- | :--- |
| `raw data` | `/BioII/chenxupeng/student/data/raw_data/*.fastq` |

**Software/Parameters:**

`cutadapt`: cutadapt removes adapter sequences from high-throughput sequencing reads.

Usage: `cutadapt -a ADAPTER [options] [-o output.fastq] input.fastq`

| `options with Parameter Setting` | function |
| :--- | :--- |
| `-u -100` | remove last 100nt so that the first 50nt is kept |
| `-q 30,30` | read quality need to be above 30 |
| `-m 16` | reads less than 15nt are removed |
| `-a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC` | cut adapt |
| `--trim-n` | trim N's on ends of reads. |

**Output:**

`*.cutadapt.fastq`

* **step three - QC after Trim**

输入文件是trim后的数据，过程与step one相同

#### **2.2.3\) Clean rRNA reads**

bowtie2可以将`.fastq`文件比对到rRNA index上从而得到**不含rRNA reads的**`.fastq`**文件以及map到rRNA index上的**`.sam`**文件**

**Input:**

2.2.2 操作结束后的`*.cutadapt.fastq`

**Software/Parameters:**

bowtie2可以Clean rRNA reads得到不含rRNA reads的`.fastq`文件以及map到rRNA index上的`.sam`文件

```text
bowtie2 -p 4 [options] -x <bt2-idx> --un <address of unmapped reads> $input_file [-S <sam>]
```

| `options with Parameter Setting` | function |
| :--- | :--- |
| `--sensitive-local(default)` | allow no mismatch, etc |
| `--norc` | do not align reverse-complement version of read |
| `--no-unal` | suppress SAM records for unaligned reads |
| `--un` `<path to unmapped reads>` | store unmapped reads |
| `-x` `<path to index>/rRNA` | indexed genome/transcriptome |
| `-S` `<path to output file>` | output file foramt as sam |

对于那些map到rRNA index上的`.sam`文件，可以用`rsem-tbam2gbam`命令转化为`.bam`文件。

```text
rsem-tbam2gbam <bt2-idx> <sam> genome_bam_output
```

**Output:**

不含rRNA reads的`.fastq`文件`*.no_rRNA.fq`，位于fastq文件夹下，详见Data Structure

map到rRNA index上的`*.<rRNA>.sam`文件，位于sam文件夹下

以及`*.<rRNA>.rsem.clean.bam`文件，位于rsem\_bam文件夹下

#### **2.2.4\) Sequential Mapping**

这步的目的就是得到比对到各种RNA类型（例如miRNA, piRNA, Y RNA和srp RNA等等）的index后得到的`.sam`文件，mapping的过程就类似于clean rRNA reads的过程。

只不过，4.2.2.3 c. 比对的index是rRNA，这里只需要 1）**把index替换成其他类型的index**，2）**将上一步比对得到的**`*.no_<some type of RNA>.fq`**作为input**，重复1）2），直至比对完所有类型至这样就可以得到各种RNA类型的比对结果。

**Input:**

`*.no_<some type of RNA>.fastq`

**Software/Parameters:**

类似2.2.3，只需修改index和input：

```text
bowtie2 -p 4 [options] -x <bt2-idx> --un <address of unmapped reads> $input_file [-S <sam>]
```

| `Parameter Setting` |
| :--- |
| `-x` `<path to index>/<some type of RNA>` |
| `*.<RNA --un by the previous step>.fq`   as `$input_file` |
| `-un<path to output>/*.no_<some type of RNA>.fq` |
| `-S` `<path to .sam file>/<some type of RNA>.sam` |

对于那些map到 index上的`.sam`文件，可以用`rsem-tbam2gbam`命令转化为`.bam`文件。

对于map到hg38上的`.sam`文件，可以用samtools的view功能转化为`.bam`文件，具体可以敲入`samtools view -h`查看怎么转化

**Output:**

不含某类型RNA reads的`.fastq`文件`*.<some type of RNA>.unAligned.fastq` \(命名方式区别于去rRNA过程，读者可根据自己习惯命名\)

map到某类型RNA index上的`*.<some type of RNA>.sam`文件

以及`*.<some type of RNA>.rsem.clean.bam`文件

**提醒:**

* reads依次map到各种类型的RNA index上，推荐次序为，`miRNA、piRNA、Y_RNA、srpRNA、tRNA、snRNA、snoRNA、lncRNA、mRNA、tucp`，最后是`hg38other`
* map的最后一步非常特殊，1）index不再是RNA\_index，是hg38,不在RNA\_index文件夹下，需要注意 2）sam转bam工具也有所不同，输出文件理论上不再应该是`*.hg38other.rsem.clean.bam`而是`*.hg38other.bam`，但是文件名的设置会影响后续代码简洁性，需要注意

**2.2.5\) length & ratio**

对mapping到不同RNA类型的index的reads，我们可以统计其长度，观察其不同RNA类型的长度分布；我们还可以统计不同RNA类型的reads所占比例，作为sample QC的参考。

**length**

这里提供统计长度的.sh，脚本位置在`/BioII/chenxupeng/student/bin/length.sh`

该脚本得到的包含长度信息的文件可以用python作格式精简处理，

```text
import pandas as pd
def get_length_table(samplename):
    '''
    sample name: Sample_N14
    '''
    pd.read_table('/BioII/chenxupeng/student/data/other_annotations/length/'+samplename+'.lengthN.stat.tsv')
    df = pd.read_table('/BioII/chenxupeng/student/data/other_annotations/length/'+samplename+'.lengthN.stat.tsv')
    df = df.pivot(index='type',columns='len',values='num')
    df = df.fillna(0)
    return df
get_length_table('Sample_N14')
```

**ratio**

这里提供统计比例的.sh脚本，位置在`/BioII/chenxupeng/student/bin/ratio.sh`

该脚本得到的包含比例信息的文件可以用python作格式精简处理，

```text
def get_counts(samplename):
    '''
    samplename: Sample_N14
    '''
    df = pd.read_table('/BioII/chenxupeng/student/data/other_annotations/counts/'+samplename+'.readsN.stat.tsv',
              names=['sample','method','type','counts']).pivot(index='type',columns='sample',values='counts')
    return df
get_counts('Sample_N14')
```

[**其他参考教程**](https://lulab.gitbook.io/training/part-ii.-basic-bioinfo-analyses/1.mapping-annotation-and-qc)

## 3\) Construct Expression Matrix 指南   

完成五个样本`Sample_N1, Sample_N7, Sample_N13, Sample_N19, Sample_N25`的expression matrix构建工作，使用mapping产生的bam文件，使用`Sample_N1, Sample_N7`的counts检查mapping和construct expression matrix是否有误。

### 3.1\) Data Structure

**inputs**

| **File format** | **Information contained in file** | **File description** | **Notes** |
| :--- | :--- | :--- | :--- |
| bam | **alignments** | Produced by mapping reads to the transcriptome. | Reads are trimmed using a proprietary version of cutAdapt. We map to transcriptome for a better sensitivity \(see details in protocol and example\). |

**outputs**

| **File format** | **Information contained in file** | **File description** | **Notes** |
| :--- | :--- | :--- | :--- |
| tsv | **gene \(ncRNA\) quantifications** | Non-normalized counts. |  |

### 3.2\) Running Scripts

#### **3.2.1\) Software/Tools**

* FeatureCounts

#### 3.2.2\) FeatureCounts

对Mapping步骤得到的不同样本不同RNA类型的`<sample>.<some type of RNA>.rsem.clean.bam`文件，进行Raw Counts的统计（无需统计hg38other），结果可输出到`.../04.counts/<sample>/<sample>.<some type of RNA>.featureCounts.counts`

**Input1:**

`.../04.counts/02.mapping/*.no_<some type of RNA>/rsem_bam/<sample>.<some type of RNA>.rsem.clean.bam`

`<annotation_file>:/BioII/chenxupeng/student/data/gtf/<some type of RNA>.gtf`

**Software usage:**

```text
featureCounts  -t exon -g transcript_id -s 1 -a <annotation_file> -o <output_file> input_file1
```

**Output:**

`<sample>.<some type of RNA>.featureCounts.counts`

### 3.3\) Merge不同RNA类型的Raw Counts

上步操作我们得到不同样本不同RNA类型的Raw Counts，现在要将这些文件合并为一个文件，代码位置在`/BioII/chenxupeng/student/bin/merge.sh`。

`proj_exRNA.featureCounts.counts.merged.mx`就是我们需要的文件

### 3.4\) 检查结果正确性

用`Sample_N1, Sample_N7`的expression matrix数据和`/BioII/chenxupeng/student/data/expression_matrix/GSE71008.txt`中相应的两个样本的参考数据计算相关系数以检查结果。可以使用pearsonr correlation coefficient衡量相关性。

$$PCC = \frac{cov(X,Y)}{\sigma X \sigma Y}$$

```python
from scipy.stats import pearsonr
pearsonr(X,Y)
```

python参考代码位于`/BioII/chenxupeng/student/bin/corr.py`

[**其他参考教程**](https://lulab.gitbook.io/training/part-ii.-basic-bioinfo-analyses/2.expression-matrix)

