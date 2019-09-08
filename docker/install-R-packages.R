# From CRAN
install.packages(c('devtools', 'VGAM', 'argparse', 'magrittr', 'readr', 'mvoutlier', 
    'ggpubr', 'fastqcr', 'ggfortify'), Ncpus=4, ask=FALSE)
# From Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks", version = "3.8")
#source('https://bioconductor.org/biocLite.R')
BiocManager::install(c('SingleCellExperiment', 'scater', 'scran', 'SCnorm',
    'EDASeq', 'RUVSeq', 'DESeq2', 'edgeR', 'sva', 'apeglm', 'gglasso', 'ggbio'), Ncpus=4, ask=FALSE, update=TRUE)
# From R-forge
install.packages('countreg', repos = c('http://download.r-forge.r-project.org',
    'https://mirrors.tuna.tsinghua.edu.cn/CRAN/'), dep = TRUE)
# From GitHub
#library(devtools)
library(remotes)
install_github('ChenMengjie/VIPER')
install_github('kassambara/easyGgplot2')
install_github("Vivianstats/scImpute")
install_github("hemberg-lab/scRNA.seq.funcs")
install_github('theislab/kBET')