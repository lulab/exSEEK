# From CRAN
install.packages(c('devtools', 'argparse', 'magrittr', 'readr', 'mvoutlier', 'ggfortify', 'SIS', 'latticeExtra'), Ncpus=4, ask=FALSE)
# From Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
#BiocManager::install("TCGAbiolinks", version = "3.8")
#source('https://bioconductor.org/biocLite.R')
BiocManager::install(c('SingleCellExperiment', 'scater', 'scran', 
    'EDASeq', 'RUVSeq', 'DESeq2', 'edgeR', 'sva'), Ncpus=4, ask=FALSE, update=TRUE)

# From GitHub
#library(devtools)
library(remotes)
#install_github('ChenMengjie/VIPER')
#install_github('kassambara/easyGgplot2')
#install_github("Vivianstats/scImpute")
install_github("hemberg-lab/scRNA.seq.funcs")
install_github('theislab/kBET')
