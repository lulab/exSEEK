## Input files
./
├── data
│   ├── annotation
│   │   └──  gtf
│   │       ├── lncRNA.gtf
│   │       ├── miRNA.gtf
│   │       ├── mRNA.gtf
│   │       └── hg38.gtf
│   ├── fastq
│   │   ├── sample_1.fastq
│   │   ├── sample_2.fastq
│   │   └── sample_3.fastq
│   ├── chrom_sizes
│   │   ├── lncRNA
│   │   ├── miRNA
│   │   ├── mRNA
│   │   └── hg38
├── metadata
│   ├── sample_classes.dataset.txt
│   └──  batch_info.dataset.txt
├── output
│   ├── tbam
│   │   └── sample_1
│   │       ├── lncRNA.bam
│   │       ├── mRNA.bam
│   │       └── miRNA.bam
│   ├── gbam
│   │   └── sample_1
│   │       ├── lncRNA.bam
│   │       ├── mRNA.bam
│   │       └── miRNA.bam
│   ├── feature_selection
│   │   └── dataset
│   │       └── Normal-HCC
│   │           └── model1
│   │               ├── features.txt
│   │               ├── metrics.txt
│   │               ├── feature_importances.txt
│   │               ├── best_model.pkl
│   │               └──  evaluation.h5
│   ├── report
│   │   ├── dataset
│   │   │   ├── dataset.html

