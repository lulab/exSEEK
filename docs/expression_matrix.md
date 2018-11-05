# Construct Expression Matrix

- full length expression matrix
- domain feature expression matrix (for small RNA)

## Data Structure

**inputs**

| **File format** | **Information contained in file** | **File description** | **Notes** |
| :--- | :--- | :--- | :--- |
| bam | **alignments** | Produced by mapping reads to the transcriptome. | Reads are trimmed using a proprietary version of cutAdapt. We map to transcriptome for a better sensitivity \(see details in protocol and example\). |

**outputs**

| **File format** | **Information contained in file** | **File description** | **Notes** |
| :--- | :--- | :--- | :--- |
| tsv | **gene \(ncRNA\) quantifications** | Non-normalized counts. |  |

## Running Scripts

### Software/Tools 

* FeatureCounts

### FeatureCounts

Make use of `<sample>.<some type of RNA>.rsem.clean.bam` file of different samples and different types of RNA as outputs of the mapping step, and the statistics of Raw Counts are performed (no need to count hg38other), with the result output to `.../04.counts/<sample>/<sample>.<some type of RNA>.featureCounts.counts`

**Input1:**

`.../04.counts/02.mapping/*.no_<some type of RNA>/rsem_bam/<sample>.<some type of RNA>.rsem.clean.bam`

`<annotation_file>:/BioII/chenxupeng/student/data/gtf/<some type of RNA>.gtf`

**Software usage:**

```text
featureCounts  -t exon -g transcript_id -s 1 -a <annotation_file> -o <output_file> input_file1
```

**Output:**

`<sample>.<some type of RNA>.featureCounts.counts`

## Merge Raw Counts of different types of RNA

In the previous step we got Raw Counts for different types of RNA from different samples. Now we need to merge these files into one file with the code location at `/BioII/chenxupeng/student/bin/merge.sh`.

`proj_exRNA.featureCounts.counts.merged.mx` is the file we need

## Check result correctness


The correlation coefficient was calculated by using the expression matrix data of `Sample_N1, Sample_N7` and the reference data of the corresponding two samples in `/BioII/chenxupeng/student/data/expression_matrix/GSE71008.txt` to check if the result is right. Correlation can be measured using the pearsonr correlation coefficient.

$$PCC = \frac{cov(X,Y)}{\sigma X \sigma Y}$$

```python
from scipy.stats import pearsonr
pearsonr(X,Y)
```
The reference code is located at `/BioII/chenxupeng/student/bin/corr.py`

## Domain Feature to Construct Expression Matrix (option)

Since the data we use is small RNA sequencing data, in addition to using full length data to construct the expression matrix, another idea is to use the peak calling method to obtain the domain feature as the expression matrix. Use the tool [piranha](https://github.com/smithlabcode/piranha) to call peak.

We provide the `Snakefile` file, input the `bam` file as outputs of mapping, and output it as the constructed expression matrix. Just ensure that the path of the input file is correct, and then execute:

```
snakemake -s Snakefile
```


