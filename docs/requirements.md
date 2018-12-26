# Requirements

## Software

* Python 3.6 (miniconda)
* Python 2.7 (miniconda)
* Java 8
* R 3.4 (https://mran.microsoft.com/download)

## Configure conda channels
```bash
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/r/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/mro/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/pro/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/
conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/
```

## Install Python packages using conda
```bash
conda install -y numpy scipy scikit-learn 
conda install -y pandas matplotlib seaborn
conda install -y tqdm snakemake h5py bokeh
conda install -y umap jinja2
```

## Install Bioconda packages

List of all available Bioconda packages: (https://bioconda.github.io/recipes.html)

```bash
conda install -y bedtools samtools star subread bowtie2
conda install -y rsem bamtools cutadapt picard gffread gffcompare
conda install -y ucsc-bedtogenepred ucsc-genepredtogtf ucsc-bedgraphtobigwig ucsc-bigwigtobedgraph
conda install -y htseq fastx_toolkit biopython rpy2
conda install -y flexbar
```

## Install Ubuntu packages

```bash
sudo apt-get install -y gzip pigz openjdk-8-jdk
```

## Install R packages

Install by running the following code in an R interactive session:
```R
options("repos" = c(CRAN="https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# From CRAN
install.packages(c('devtools', 'sva', 'VGAM', 'argparse', 'magrittr', 'readr', 'mvoutlier', 'ggpubr'))
# From Bioconductor
source('https://bioconductor.org/biocLite.R')
biocLite(c('SingleCellExperiment', 'scater', 'scran', 'SCnorm',
    'EDASeq', 'RUVSeq', 'DESeq2', 'edgeR', 'sva', 'apeglm'))
# From R-forge
install.packages('countreg', repos = c('http://R-Forge.R-project.org', 'https://mirrors.tuna.tsinghua.edu.cn/CRAN/'), dep = TRUE)
# From GitHub
library(devtools)
install_github('ChenMengjie/VIPER')
install_github('kassambara/easyGgplot2')
install_github("Vivianstats/scImpute")
install_github("hemberg-lab/scRNA.seq.funcs")
```

## Other packages
* [find_circ 1.2](https://github.com/marvin-jens/find_circ) (depends on Python 2.7)

## Singularity

### Build image

```bash
singularity build singularity/exseek.img singularity/Singularity
```

### Make wrappers for singularity executables
```bash
bin/make_singularity_wrappers.py \
    --image ~/singularity/simg/exseek.simg \
    --list-file singularity/exports.txt \
    --singularity-path $(which singularity) \
    -o ~/singularity/wrappers/exseek
```

### Add wrappers to PATH
```bash
export PATH="$HOME/singularity/wrappers/exseek:$PATH"
```
