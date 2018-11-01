
install.packages(c('SingleCellExperiment', 'scater', 'scImpute', 'scran', 'SCnorm', 'EDASeq', 'RUVSeq', 'sva', 'scRNA.seq.funcs'))

BiocInstaller::biocLite(c('SingleCellExperiment', 'scater', 'scImpute', 'scran', 'SCnorm', 'EDASeq', 'RUVSeq', 'sva', 'scRNA.seq.funcs'))


lapply(c('SingleCellExperiment', 'scater', 'scImpute', 'scran', 'SCnorm', 'EDASeq', 'RUVSeq', 'sva', 'scRNA.seq.funcs'), function(x){if (!(x %in% .packages(T))) BiocInstaller::biocLite(x)})

# edit .Rprofile
install.packages('R.oo')
lapply(c('SingleCellExperiment', 'scater', 'scImpute', 'scran', 'SCnorm', 'EDASeq', 'RUVSeq', 'sva', 'scRNA.seq.funcs'), function(x){if (!(x %in% .packages(T))) BiocInstaller::biocLite(x)})
