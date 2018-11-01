# edit .Rprofile
install.packages('R.oo')
lapply(c('SingleCellExperiment', 'scater', 'scran', 'SCnorm', 'EDASeq', 'RUVSeq', 'sva'), function(x){if (!(x %in% .packages(T))) BiocInstaller::biocLite(x)})

devtools::install_github('hemberg-lab/scRNA.seq.funcs')
devtools::install_github('Vivianstats/scImpute')
