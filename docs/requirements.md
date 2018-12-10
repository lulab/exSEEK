**Software**

* Python 3.6 (miniconda)
* Python 2.7 (miniconda)
* Java 8
* R 3.4

**Install Python packages using conda**
```bash
conda install -y numpy scipy scikit-learn pandas matplotlib seaborn tqdm snakemake h5py bokeh
```

**Install Bioconda packages**

List of all available Bioconda packages: (https://bioconda.github.io/recipes.html)

```bash
conda install -y bedtools samtools STAR featureCounts bowtie2 \
    rsem bamtools cutadapt picard gffread gffcompare \
    ucsc-bedtogenepred ucsc-genepredtogtf ucsc-bedgraphtobigwig \
    htseq fastx_toolkit biopython rpy2
```

**Install Ubuntu packages**

```bash
sudo apt-get install -y gzip pigz
```

**Install R packages**

Install by running the following code in an R interactive session:
```R
# From CRAN
install.packages(c('devtools', 'SingleCellExperiment', 'scater', 'scImpute', 'scran', 'SCnorm',
    'EDASeq', 'RUVSeq', 'sva', 'scRNA.seq.funcs', 'VGAM', 'argparse', 'magrittr'))
install.packages(c('devtools', 'sva', 'VGAM', 'argparse', 'magrittr'))
# From Bioconductor
source('https://bioconductor.org/biocLite.R')
biocLite(c('SingleCellExperiment', 'scater', 'scImpute', 'scran', 'SCnorm',
    'EDASeq', 'RUVSeq', 'scRNA.seq.funcs'))
# From R-forge
install.packages('countreg', repos = c('http://R-Forge.R-project.org', 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/'), dep = TRUE)
# From GitHub
library(devtools)
install_github('ChenMengjie/VIPER')
install_github('kassambara/easyGgplot2)
```

**Other packages**
* [find_circ 1.2](https://github.com/marvin-jens/find_circ) (depends on Python 2.7)