#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()
parser$add_argument("-s", "--step", required=TRUE, help="which step to run")
parser$add_argument("-i", "--input", required=TRUE, help="input expression matrix file")
parser$add_argument("-o", "--output", required=TRUE, help="output expression matrix file")
parser$add_argument("--temp-dir", required=FALSE, default=".", help="temporary directory")
parser$add_argument("-c", "--class", required=FALSE, help="input class info file")
parser$add_argument("-b", "--batch", required=FALSE, help="input batch info file")
parser$add_argument("-m", "--method", required=FALSE, help="name of the method")

parser$add_argument("--filtercount", type="integer", default=10, help="count threshold for filtering")
parser$add_argument("--filtercpm", type="double", default=10, help="CPM threshold for filtering")
parser$add_argument("--filterrpkm", type="double", help="RPKM threshold for filtering")

parser$add_argument( "--filtermethod", type="character", default="filtercount",
                    metavar="STRING",
                    help="the filter algorithm to use [default = %(default)s], and return count matrix")
parser$add_argument("--filterexpv", type="double", default=5,
    help="filter by expression value of a gene [default = %(default)s]",
    metavar="NUMBER")

parser$add_argument("--filtersample", type="double", default=0.2,
    help="filter by counts of sample above certain counts of a gene [default = %(default)s]",
    metavar="NUMBER")
parser$add_argument("--imputecluster", type="integer", default=5,
    help="cluster number in scImpute [default = %(default)s]",
    metavar="NUMBER")
parser$add_argument("--imputevipernum", type="integer", default=5000,
    help="number in viper [default = %(default)s]",
    metavar="NUMBER")
parser$add_argument( "--imputecutoff", type="double", default=0.5,
                    metavar="NUMBER",
                    help="cutoff in viper [default = %(default)s]")
parser$add_argument( "--imputealpha", type="double", default=0.1,
                    metavar="NUMBER",
                    help="alpha in viper [default = %(default)s]")

parser$add_argument( "--normtopk", type="integer", default=20,
                    metavar="NUMBER",
                    help="top K feature as scale factor [default = %(default)s]")
parser$add_argument( "--cvthreshold", type="double", default=0.5,
                    metavar="NUMBER",
                    help="coefficient variance threshold of reference gene, filter ref gene with CV bigger than [default = %(default)s]")
parser$add_argument( "--remove-gene-types", type="character", default="miRNA,piRNA",
                    metavar="STRING",
                    help="remove some time of RNA for normalization scale factor calculation [default = %(default)s]")
parser$add_argument( "--ref-gene-file", type="character", #default="miRNA,piRNA",
                    metavar="STRING",
                    help="reference gene file path [default = %(default)s]")
#they are feature name for full length feature, most are miRNA, for domain feature, they have the same feature name

parser$add_argument("--batch-index", type="integer", default=1,
                    metavar="INT",
                    help="batch index to select which batch to use [default = %(default)s]")
parser$add_argument("--ruv-k", type="integer", default=1, metavar="INT", 
    help="parameter k for RUVs (The number of factors of unwanted variation to be estimated from the data)")

parser$add_argument("-p", "--processors", type="integer", default=1,
    help="Number of processors to use. This option is useful on multicore *nix or Mac machine only, when performing multiple runs (nrun > 1) [default = %(default)s]",
    metavar="NUMBER")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults,
args <- parser$parse_args()

splitname <- unlist(strsplit(args$input,split = '/',fixed = TRUE))
splitname <- splitname[length(splitname)]

#' @title read counts matrix
#'
#' @param path string.
#' @param ... other arguments passsed on to [readr::read_tsv()]
#'
#' @return integer matrix
#'
#' @details In any case, first part (separated by `|`) of row names must be
#'   Ensembl transcript id
#'
#' @export


#' @title sample classinfo
#'
#' @param path string.
#'
#' @return string matrix
#'
#' @details column 1 represents sample name, column 2 represents classinfo
#'
#' @export

# path = 'scirep_classes.txt'
read_classinfo <- function(path, ...) {
    read.table(path, sep='\t', header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
}

read_matrix <- function(filename){
    read.table(filename, sep='\t', header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
}

write_matrix <- function(mat, filename){
    write.table(mat, filename, sep='\t',quote=FALSE, row.names=TRUE, col.names=TRUE)
}

################################################################################
#################################imputation#####################################
################################################################################

#' @title filter genes with low expression values
#'
#' @param mat integer matrix.
#' @param min_count, min_sample_per_gene integer scalar. For each gene, it must
#'   contain at least `min_count` reads in at least `min_sample_per_gene`
#'   samples. Otherwise, it would be dropped.
#'
#' @return integer matrix.
#'
#' @examples
#' filter_low(sim_mat)
#'
#' @export
#filter_low <- function(mat, min_count = 5, min_sample_per_gene = 10) {
#    print(paste('start filtering lowly expressed gene:','count threshold',min_count,'sample threshold',min_sample_per_gene,sep=' '))
#	low_per_row <- rowSums(mat > min_count)
#	keeped_row <- low_per_row > min_sample_per_gene
#	mat[keeped_row, ]
#}

filter_low <- function(mat, min_count = 5, min_sample_per_gene = 0.5) {
    print(paste('start filtering lowly expressed gene:','count threshold',min_count,'sample threshold',min_sample_per_gene,sep=' '))
    min_sample_per_gene <- ceiling(dim(mat)[2]*min_sample_per_gene)
    low_per_row <- rowSums(mat > min_count)
    keeped_row <- low_per_row >= min_sample_per_gene
    mat[keeped_row, ]
}

filter_low_cpm <- function(mat, min_cpm = 5, min_pct_sample_per_gene = 0.5) {
    print(paste('start filtering lowly expressed gene:','CPM threshold',min_cpm,'sample threshold',min_pct_sample_per_gene,sep=' '))
    row_all <- nrow(mat) %>% seq_len()
    mat_cpm <- t(t(mat*1e6) / colSums(mat[row_all, , drop = F], na.rm = T))
    min_sample_per_gene <- ceiling(dim(mat_cpm)[2]*min_pct_sample_per_gene)
    low_per_row <- rowSums(mat_cpm > min_cpm)
    keeped_row <- low_per_row >= min_sample_per_gene
    mat[keeped_row, ]
}

filter_low_rpkm <- function(mat, min_rpkm = 5, min_pct_sample_per_gene = 0.5) {
    print(paste('start filtering lowly expressed gene:','RPKM threshold',min_rpkm,'sample threshold',min_pct_sample_per_gene,sep=' '))
    row_all <- nrow(mat) %>% seq_len()
    mat_cpm <- t(t(mat*1e6) / colSums(mat[row_all, , drop = F], na.rm = T))
    
    gene_length <- c()
    for(i in seq_len(length(rownames(mat_cpm)))){
            gene_length[i] <- as.integer(unlist(strsplit(rownames(mat_cpm)[i],"|",fixed=T))[7])
                             -as.integer(unlist(strsplit(rownames(mat_cpm)[i],"|",fixed=T))[6])
    }
    mat_rpkm <- mat_cpm*1000/gene_length
    
    min_sample_per_gene <- ceiling(dim(mat_rpkm)[2]*min_pct_sample_per_gene)
    low_per_row <- rowSums(mat_rpkm > min_rpkm)
    keeped_row <- low_per_row >= min_sample_per_gene
    mat[keeped_row, ]
}

#' @imputation
#'
#' @param mat integer matrix.
#' @param tmp_path where tmp files stores, "data/expression_matrix/" for example.
#' @param out_path where outputs stores, "data/matrix_processing/imputation/" for example.
#' @param K imputation Kcluster
#' @param N imputation ncores
#' @return integer matrix named "scimpute_count.txt" stored in out_path.
#'
#' @examples imputation(mat, "data/expression_matrix/", "data/matrix_processing/imputation/",5,3)
#'
#' @export
scimpute_count <- function(mat, temp_dir, K = 5, N = 3) {
    suppressMessages(library("scImpute"))
    print('start imputation using scImpute')
    mat_correct <- names(mat)
    names(mat) <-  paste('C_',seq_len(length(names(mat))))
    write.csv(mat, paste(temp_path, "input.csv",sep=""), sep=',')
    dir.create(splitname)
    scimpute(count_path = paste(temp_dir, "input.csv",sep=""), infile = "csv", outfile = "txt", out_dir = temp_dir, Kcluster = K, ncores = N)
    mat <- read.table(paste(temp_dir, "scimpute_count.txt", sep="/"),sep=' ',header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
    unlink(temp_dir, recursive=TRUE)
    names(mat) <-mat_correct
    mat
}

viper_count <- function(mat, num = 5000, percentage.cutoff = 0.1, alpha= 0.5, temp_dir=".") {
    suppressWarnings(library(VIPER))
    mat_correct <- names(mat)
    names(mat) <-  paste('C_',seq_len(length(names(mat))))
    print (num, percentage.cutoff, alpha)
    VIPER(mat, num = num, percentage.cutoff = percentage.cutoff, minbool = FALSE, alpha = alpha, report = TRUE, outdir = temp_dir)
    mat <- read.table(paste(temp_dir, 'imputed_counts.csv', sep='/'),sep=' ',header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
    unlink(temp_dir, recursive=TRUE)
    names(mat) <-mat_correct
    mat
}


################################################################################
###############################normalization####################################
################################################################################
#suppressPackageStartupMessages(library(clusterSim))
suppressPackageStartupMessages(library(scRNA.seq.funcs))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(kBET))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(ggpubr))
options(stringsAsFactors = FALSE)
library(magrittr)

#' @title martix normalization
#' @examples
#' \donotrun{
#'     norm_mat(
#'         '/path/to/matrix'
#'     )
#' }
normalize <- function(
    mat,
    method,
    top_n = 20, 
    rm_gene_types = NULL,
    ref_genes=NULL
) {
    if (method == 'SCnorm')    mat <- norm_SCnorm(mat)
    else if (method == 'TMM')       mat <- norm_tmm(mat)
    else if (method == 'RLE')       mat <- norm_rle(mat)
    else if (method == 'CPM')       mat <- norm_cpm_total(mat)
    else if (method == 'UQ')        mat <- norm_uq(mat)
    else if (method == 'CPM_top')   mat <- norm_cpm_top(mat, top_n)
    else if (method == 'CPM_rm')    {
        if(is.null(rm_gene_types)){
            stop('argument rm_gene_types is required for normalization method: CPM_top')
        }
        mat <- norm_cpm_rm(mat, rm_gene_types)
    }
    else if (method == 'CPM_refer') mat <- norm_cpm_refer(mat, refer_gene_id_path)
    else if (method == 'null')      mat <- mat
    else stop('unknown normalization method: ', method)
    mat
}


#' @title SCnorm normalization
#'
#' @param mat integer matrix. counts
#' @param ... other arguments passed on to [SCnorm::SCnorm()]
#'
#' @examples
#' norm_SCnorm(sim_mat*10)
#'
#' @family matrix normalization
#'
#' @export
norm_SCnorm <- function(mat, ...) {
    print('start normalization using SCnorm')
    Conditions = rep(1, ncol(mat));
    sce <- suppressMessages(SCnorm::SCnorm(mat, Conditions, NCores=4));
    SCnorm::results(sce)
}

# norm_scater ------------------

#' @title TMM/RLE normalization by scater package
#'
#' @param mat integer matrix. counts
#'
#' @family matrix normalization
#'
#' @name norm_scater

as_SingleCellExperiment <- function(mat, col_data = NULL) {
    assays = list(counts = as.matrix(mat))
    if (is.null(col_data))
    SingleCellExperiment::SingleCellExperiment(assays = assays)
    else
    SingleCellExperiment::SingleCellExperiment(assays = assays, colData = col_data)
}

#' @rdname  norm_scater
#'
#' @details `norm_uq()` performs upper-quartile normalization
#'
#' @examples
#' norm_uq(sim_mat)
#'
#' @export
norm_uq <- function(mat) {
    print('start normalization using UQ')
    #mat %>% as_SingleCellExperiment() %>%
    #{suppressWarnings(scater::normaliseExprs(., "TMM"))} %>%
    #scater::normalise() %>% SingleCellExperiment::normcounts()
    dl <- edgeR::DGEList(counts=mat)
    dl <- edgeR::calcNormFactors(dl, method='upperquartile')
    edgeR::cpm(dl)
}

#' @rdname  norm_scater
#'
#' @details `norm_tmm()` performs TMM normalization
#'
#' @examples
#' norm_tmm(sim_mat)
#'
#' @export
#norm_tmm <- function(mat) {
#    print('start normalization using TMM')
#    mat %>% as_SingleCellExperiment() %>%
#    {suppressWarnings(scater::normaliseExprs(., "TMM"))} %>%
#    scater::normalise() %>% SingleCellExperiment::normcounts()
#}

norm_tmm <- function(mat) {
    print('start normalization using TMM')
    dl <- edgeR::DGEList(counts=mat)
    dl <- edgeR::calcNormFactors(dl, method='TMM')
    return(edgeR::cpm(dl))
}


norm_UQ <- function(mat) {
    print('start normalization using upperquartile')
    dl <- edgeR::DGEList(counts=mat)
    dl <- edgeR::calcNormFactors(dl, method='upperquartile')
    return(edgeR::cpm(dl))
}
#' @rdname  norm_scater
#'
#' @details `norm_rle()` performs RLE normalization
#'
#' @examples
#' norm_rle(sim_mat)
#'
#' @export
#norm_rle <- function(mat) {
#    print('start normalization using RLE')
#    mat %>% as_SingleCellExperiment() %>%
#    {suppressWarnings(scater::normaliseExprs(., "RLE"))} %>%
#    scater::normalise() %>% SingleCellExperiment::normcounts()
#}
norm_rle <- function(mat) {
    print('start normalization using RLE')
    dl <- edgeR::DGEList(counts=mat)
    dl <- edgeR::calcNormFactors(dl, method='RLE')
    edgeR::cpm(dl)
}

# norm_cpm ------------------

#' @title CPM normalization by some genes
#'
#' @param mat integer matrix. counts
#' @param row integer or logical. Use which rows (genes) as normalization factor
norm_cpm_impl <- function(mat, row) {
    t(t(mat*1e6) / colSums(mat[row, , drop = F], na.rm = T))
    #edgeR::cpm(mat)
}


#' @title CPM normalization
#'
#' @description CPM normalization using counts sum of _certain_ genes as scaling factor
#'
#' @param mat integer matrix. counts.
#'
#' @details some functions may throw errors
#'
#' @family matrix normalization
#'
#' @name norm_cpm




#' @rdname norm_cpm
#'
#' @details `norm_cpm_total()` uses total genes
#'
#' @examples
#' norm_cpm_total(sim_mat)
#'
#' @export
norm_cpm_total <- function(mat) {
    print('start normalization using CPM')
    row_all <- nrow(mat) %>% seq_len()
    
    norm_cpm_impl(mat, row_all)
}

#' @rdname norm_cpm
#'
#' @param top_n integer scalar. see `norm_cpm_top()` below
#'
#' @details `norm_cpm_top()` uses top 20 genes sorted by counts (assuming `top_n = 20L`)
#'
#' @examples
#' norm_cpm_top(sim_mat, 20L)
#'
#' @export
#norm_cpm_top <- function(mat, top_n) {
#   print(paste('start normalization using top',top_n,'genes as scale factor',sep=' '))
#   if (nrow(mat) < top_n)
#    stop('too few feature for CPM top n normalization')
#    
#    row_top <-  mat %>% rowSums() %>% sort(decreasing = T, index.return = T) %>%
#    {.$ix[seq_len(top_n)]}
#    
#    norm_cpm_impl(mat, -row_top)
#}
norm_cpm_top <- function(mat, top_n) {
    print(paste('start normalization using top',top_n,'genes as scale factor',sep=' '))
    if (nrow(mat) < top_n)
    stop('too few feature for CPM top n normalization')
    
    row_top <-  mat %>% rowSums() %>% sort(decreasing = T, index.return = T) %>%
    {.$ix[seq_len(top_n)]}
    
    top = t(t(mat[row_top,]*1e6) / colSums(mat[row_top, , drop = F], na.rm = T))
    top_down= t(t(mat[setdiff(seq_len(dim(mat)[1]),row_top),]*1e6) / colSums(mat[setdiff(seq_len(dim(mat)[1]),row_top), , drop = F], na.rm = T))
    mat_top <- rbind(top,top_down)
    mat_top[rownames(mat),]
}

#' @rdname norm_cpm
#'
#' @param gene_type character. see `norm_cpm_rm()` below
#'
#' @details `norm_cpm_rm()` uses non-piRNA genes (assuming `gene_type = 'piRNA'`)
#'
#' @examples
#' norm_cpm_rm(sim_mat, c('miRNA', 'piRNA'))
#'
#' @export
norm_cpm_rm <- function(mat, rm_gene_type) {
        print(paste('start normalization by removed some kind of RNA type ',args$removetype,sep=' '))
        row_rm <-  mat %>% rownames() %>% strsplit(split='|',fixed=TRUE) %>% data.frame() %>% {.[2,]} %>% {. %in% rm_gene_type}
        return(norm_cpm_impl(mat, !row_rm))
}



#' @rdname norm_cpm
#'
#' @param refer_gene_id character. Ensembl transcript id, see `norm_cpm_refer()` below
#'
#' @details `norm_cpm_refer()` uses given reference genes
#'
#' @examples
#' norm_cpm_refer(sim_mat, suggest_refer$id)
#'
#' @export

# mat = sim_mat
# refer_gene_id = suggest_refer$id
cv_fun <- function(x) {
	sd(x, na.rm = T) / mean(x, na.rm = T)
}


norm_cpm_refer <- function(mat, ref_genes, cv_threshold=0.5) {
    message('start normalization by reference gene with CV threshold', cv_threshold)
    keeped_ref <- mat[ref_genes, , drop = F] %>% apply(1, cv_fun) < cv_threshold
    if(!any(kept_ref)){
        stop('no reference genes left after filtering by CV')
    }
    norm_cpm_impl(mat, ref_genes[keeped_ref])
}

################################################################################
#################################batch removal##################################
################################################################################


remove_batch <- function( mat, method, class_info=NULL, batch_info=NULL, ruv_k=1){
    # only remove batch for samples with batch information
    if(!is.null(batch_info)){
        samples_with_batch <- names(batch_info)[!is.na(batch_info)]
    }else{
        samples_with_batch <- colnames(mat)
    }
    if (method == 'RUV')       mat <- ruvs(mat, class_info=class_info, k=ruv_k)
    else if(method == 'RUVn')       {
        class_info <- rep(1, ncol(mat))
        names(class_info) <- colnames(mat)
        mat <- ruvs(mat, class_info=class_info, k=ruv_k)
    }
    else if (method == 'ComBat')    mat[, samples_with_batch] <- combat(mat[, samples_with_batch], class_info=class_info, batch_info=batch_info)
    else if (method == 'limma')     mat[, samples_with_batch] <- limma(mat[, samples_with_batch], class_info=class_info, batch_info=batch_info)
    else if (method == 'null')      mat <- mat
    else stop("unknown batch effect removal method: ", method)
    mat
}

ruv <- function(
    mat,
    classinfo_path,
    label_column = 2,
    k = 10
){
    suppressMessages(library(RUVSeq))

    print('start batch removal using RUVs')
    cIdx <- rownames(mat)
    
    sample_info <- read.table(classinfo_path,sep='\t',header=TRUE,  check.names=FALSE,  stringsAsFactors=FALSE)
    ##rank by mat
    if(unique(is.na(sample_info$sample_id))) 
    stop("sample_id not in file")
    rownames(sample_info) = sample_info$sample_id
    sample_info=sample_info[names(mat),]
    rownames(sample_info) <- c()
    
    names(sample_info)[label_column]="label"
    scIdx <- matrix(-1, ncol = max(table(sample_info$label)), nrow = dim(table(sample_info$label)))
    labellist <- names(table(sample_info$label))
    for(i in c(1:dim(table(sample_info$label)))) {
        tmp <- which(sample_info$label == labellist[i])
        scIdx[i, 1:length(tmp)] <- tmp
    }
    mat <- log(mat+0.001)
    ruv <- RUVs(as.matrix(mat), cIdx, k = k, scIdx = scIdx, isLog = TRUE)
    exp(ruv$normalizedCounts)
}

ruvs <- function(mat, class_info, batch_info=NULL, k = 1){
    if(is.null(class_info)) stop('class_info is needed for RUVs')

    message('start batch removal using RUVs')
    suppressMessages(library(RUVSeq))

    cIdx <- rownames(mat)
    
    class_sizes <- table(class_info)
    scIdx <- matrix(-1, ncol = max(class_sizes), nrow = dim(class_sizes))
    for(i in c(1:dim(class_sizes))) {
        tmp <- which(class_info == names(class_sizes)[i])
        scIdx[i, 1:length(tmp)] <- tmp
    }
    mat <- log(mat + 0.25)
    seq_ruvs <- RUVs(as.matrix(mat), cIdx, k = k, scIdx = scIdx, isLog = TRUE)
    exp(seq_ruvs$normalizedCounts)
}


combat <- function(mat,class_info=NULL, batch_info=NULL){
    if(is.null(batch_info)) stop('batch_info is needed for ComBat')

    message('start batch removal using combat')
    suppressMessages(library(sva))
    #batch_info <-read.table(batchinfo_path,sep='\t',row.names=1,header=T,check.names = FALSE)
    #if (!(dim(mat)[2]==dim(batch_info)[1]))
    #    stop('sample numbers in batch info and expression matrix should be same')
    #batchname <-toString(names(batch_info)[batch_column])
    #batch_info=as.data.frame(batch_info[names(mat),])
    #batch_info <- as.factor(batch_info)
    mod <- model.matrix(~ 1, data = as.factor(batch_info))
    combat <- ComBat(
        dat = log(as.matrix(mat) + 0.25),
        batch = as.factor(batch_info),
        mod = mod,
        par.prior = TRUE,
        prior.plots = FALSE
    )
    
    mat <- exp(combat)
    mat
}

limma <- function(
    mat,
    class_info=NULL,
    batch_info=NULL
){
    if(is.null(batch_info)) stop('batch_info is needed for limma')

    print('start batch removal using limma')
    suppressMessages(library(limma))

    mat <- removeBatchEffect(log(as.matrix(mat) + 0.25), as.factor(batch_info))
    mat <- exp(mat)
    mat
}

################################################################################
#################################plot part######################################
################################################################################

#' @export
plot_highest_exprs <- function(sce, top_n = 20) {
	sce %>% {suppressMessages(scater::calculateQCMetrics(.))} %>%
		scater::plotHighestExprs(n = top_n)
}


# plot_group --------------

plot_group_impl <- function(sce, shape = NULL, color = NULL, plot_fun) {
 	plot_fun(
 		sce,
		shape_by = shape, colour_by = color,
    	run_args = list(exprs_values = 'counts')
	)
}

#' @title plot PCA, TSNE
#'
#' @param sce A SingleCellExperiment object.
#' @param shape, color string. specify a column in `col_data` of [as_SingleCellExperiment()] to shape/color by
#'
#' @name plot_group



#' @rdname plot_group
#'
#' @examples
#' as_SingleCellExperiment(sim_mat) %>% plot_PCA()
#' 
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA()
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(color = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label', color = 'label')
#'
#' @export
plot_PCA <- function(sce, shape = NULL, color = NULL) {
	plot_group_impl(sce, shape, color, scater::plotPCA)
}



#' @rdname plot_group
#'
#' @examples
#' as_SingleCellExperiment(sim_mat) %>% plot_PCA()
#' 
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA()
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(color = 'label')
#' as_SingleCellExperiment(sim_mat, sim_sample_class) %>% plot_PCA(shape = 'label', color = 'label')
#'
#' @export
plot_TSNE <- function(sce, shape = NULL, color = NULL) {
	plot_group_impl(sce, shape, color, scater::plotTSNE)
}

# plot_variance --------------


#' get y axis range of ggplot object.
get_y_range <- function(plot) {
    ggplot2::ggplot_build(plot)$layout$panel_params[[1]]$y.range
}

#' @title generate equally spaced y coordinates, not hit bottom nor top.
#'
#' @details Image `plot`'s y axis extends from 0 to 3, `x` contains 3 values,
#' then we give `c(0.5, 1.5, 2.5)`.
#'
#' @param plot ggplot object.
#' @param x numberic vector
#'
#' @return numeric vector. The same length as `x`
#'
#' @keywords internal
seq_y <- function(plot, x) {
    y_range <- get_y_range(plot)
    by <- diff(y_range) / length(x)
    seq(y_range[1] + by, y_range[2], by) - by/2
}


#' @title plot variance of counts across samples
#'
#' @param mat integer matrix. counts
#' @param refer_gene_id character. Ensembl transcript id, add a vertical line
#'   for each gene to mark the corresponding CV (on x axis). Only genes in
#'   counts matrix would be shown. Usually these genes should be the reference
#'   genes you want to use for normalization.
#' @param refer_gene_name character. Transcript name
#'
#' @return [ggplot2::ggplot()] object.
#'
#' @name plot_variance






#' @details
#' `plot_cv_density()` produces density plot of coefficient of variation
#'
#' @examples
#' # only one gene exist in the matrix
#' plot_cv_density(sim_mat, suggest_refer$id)
#' plot_cv_density(sim_mat, suggest_refer$id, suggest_refer$name)
#'
#' # the name should be the same length as id
#' plot_cv_density(sim_mat, rownames(sim_mat)[1:6], letters[1:3])
#' # if only part of the genes have name, you can pass the id of other genes
#' plot_cv_density(sim_mat, rownames(sim_mat)[1:6], c(letters[1:3], rownames(sim_mat)[4:6]))
#'
#' @export
#'
#' @rdname plot_variance

# mat = sim_mat
# refer_gene_id = suggest_refer$id
# refer_gene_name = suggest_refer$name
plot_cv_density <- function(mat, refer_gene_id = '', refer_gene_name = refer_gene_id) {
    cv <- mat %>% apply(1, cv_fun) %>%
    {tibble::tibble(id = names(.), value = .)} %>%
    dplyr::mutate(id = stringr::str_extract(id, '[^|]+'))
    plot <- ggplot2::ggplot(cv, ggplot2::aes(value)) +
    ggplot2::geom_density(color = 'blue') +
    ggplot2::labs(x = 'coefficient of variation')
    
    if (length(refer_gene_id) != length(refer_gene_name)) {
        warning("Ignoring refer_gene_name, since it isn't the same length as refer_gene_id")
        refer_gene_name = refer_gene_id
    }
    cv_refer <- tibble::tibble(id = refer_gene_id, name = refer_gene_name) %>%
    dplyr::inner_join(cv, by = 'id')
    if (nrow(cv_refer) == 0L) {
        warning("None refer gene found in the count matrix")
        return(plot)
    }
    
    plot + ggplot2::geom_vline(xintercept = cv_refer$value, color = 'green') +
    ggplot2::geom_point(
    ggplot2::aes(x = value, y = seq_y(plot, value)),
    data = cv_refer, size = 2, shape = 1
    ) +
    ggrepel::geom_label_repel(
    ggplot2::aes(x = value, y = seq_y(plot, value), label = name),
    data = cv_refer, hjust = 0.5
    )
}



#' @details
#' `plot_cv_density()` produces density plot of coefficient of variation
#'
#' @examples
#' # only one gene exist in the matrix
#' plot_refer_violin(sim_mat, suggest_refer$id)
#' plot_refer_violin(sim_mat, suggest_refer$id, suggest_refer$name)
#'
#' # the name should be the same length as id
#' plot_refer_violin(sim_mat, rownames(sim_mat)[1:6], letters[1:3])
#' # if only part of the genes have name, you can pass the id of other genes
#' plot_refer_violin(sim_mat, rownames(sim_mat)[1:6], c(letters[1:3], rownames(sim_mat)[4:6]))
#'
#' @export
#'
#' @rdname plot_variance

# mat = sim_mat
# refer_gene_id = rownames(mat)[1:6]
# refer_gene_name = paste0('gene_', letters[1:6])
plot_refer_violin <- function(mat, refer_gene_id, refer_gene_name = refer_gene_id) {
    if (length(refer_gene_id) != length(refer_gene_name)) {
        warning("Ignoring refer_gene_name, since it isn't the same length as refer_gene_id")
        refer_gene_name = refer_gene_id
    }
    
    refer_gene <- tibble::tibble(id = refer_gene_id, name = refer_gene_name)
    refer_count <- mat %>% tibble::as_tibble(rownames = 'id') %>%
    dplyr::mutate(id = stringr::str_extract(id, '[^|]+')) %>%
    dplyr::inner_join(refer_gene, ., by = 'id') %>% dplyr::select(-id)
    if (nrow(refer_count) == 0L) {
        warning('None refer gene found in the count matrix')
        return(ggplot2::ggplot())
    }
    
    refer_count_long <- refer_count %>%    tidyr::gather('sample', 'count', -1) %>%
    dplyr::mutate_at('name', as.factor)
    g_violin <- refer_count_long %>%
    ggplot2::ggplot(ggplot2::aes(name, log2(count + 0.001))) +
    ggplot2::geom_violin() +
    ggplot2::labs(x = 'reference transcripts', y = quote(log[2](count)))
    
    # max y coordinate of each violin
    y_max <- ggplot2::ggplot_build(g_violin)$data[[1]] %>% tibble::as_tibble() %>%
    dplyr::group_by(x) %>% dplyr::arrange(dplyr::desc(y)) %>% dplyr::slice(1) %>%
    dplyr::ungroup() %>% dplyr::arrange(x) %>% dplyr::select(x, y)
    
    cv_df <- refer_count_long %>%
    dplyr::group_by(name) %>% dplyr::summarise(cv = cv_fun(count)) %>%
    dplyr::arrange(name) %>% dplyr::mutate(x = seq_along(name)) %>%
    dplyr::inner_join(y_max, by = 'x') %>%
    dplyr::mutate(y = y + diff(get_y_range(g_violin)) / 20) %>%
    dplyr::mutate(cv = formatC(cv, digits = 3, format = 'f'))
    
    g_violin + ggplot2::geom_text(ggplot2::aes(x, y, label = cv), cv_df, color = 'blue')
}

################################################################################
#########################process pipeline#######################################
################################################################################

#args$input: matrix_path = 'output/scirep/count_matrix/transcript.txt'
#args$class: classinfo_path = 'data/labels/scirep_classes.txt'
#args$batch: 'data/other_annotations/scirep_batch.txt'
#args$imputeout: impute_path = "output/matrix_processing/imputation/"

#dummy_io<- function(filename){}


message('read expression matrix: ', args$input)
mat <- read_matrix(args$input)
sample_ids <- colnames(mat)

if((args$step != 'filter') && (is.null(args$method))){
    stop('--method is required for step: ', args$step)
}
# filter
if (args$step =='filter'){
    if(!is.null(args$filtercount)){
        #message(sprintf('Filter features with count <= %d in %f %% samples', args$filtercount, args$filtersample*100))
        mat <- filter_low(mat, args$filtercount, args$filtersample)
    }
    else if(!is.null(args$filtercpm)){
        #message(sprintf('Filter features with CPM <= %d in %f %% samples', args$filtercpm, args$filtersample*100))
        mat <- filter_low_cpm(mat, args$filtercpm, args$filtersample)
    }
    else if(!is.null(args$filterrpkm)){
        #message(sprintf('Filter features with RPKM <= %d in %f %% samples', args$filterrpkm, args$filtersample*100))
        mat <- filter_low_rpkm(mat, args$filterrpkm, args$filtersample)
    }
} else if(args$step =='imputation'){
    # imputation
    library(readr)
    mat <-read.table(args$input,sep='\t',header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
    if (args$method == 'scimpute_count'){
        mat <- scimpute_count(mat, impute_path= args$imputeout, K = args$imputecluster, N = args$processors)
    } else if(args$method == 'viper_count'){
        mat <- viper_count(mat, impute_path= args$imputeout,num = args$imputevipernum,percentage.cutoff = args$imputecutoff, alpha= args$imputealpha) 
    } else if(args$method == 'null'){
        
    } else {
        stop('unknown imputation method: ', args$method)
    }
}  else if(args$step =='normalization'){
    # normalization
    if(!is.null(args$ref_gene_file)){
        ref_genes <- read.table(args$ref_gene_file)[,1]
    }else{
        ref_genes <- NULL
    }
    mat <- normalize(mat, method=args$method,top_n = args$normtopk, 
        rm_gene_types = args$remove_gene_types, ref_genes=ref_genes)
}  else if(args$step =='batch_removal'){
    class_info <- NULL
    batch <- NULL
    if(args$method == 'RUV'){
        if(!is.null(args$class)){
            message('read class information: ', args$class)
            class_info <- read.table(args$class, sep='\t', header=TRUE, check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
            class_info <- class_info[sample_ids, 1]
        }else{
            stop('argument --class is required for RUV')
        }
    }else if(args$method %in% c('ComBat', 'limma')){
        if(!is.null(args$batch)){
            message('read batch information: ', args$batch)
            batch <- read.table(args$batch, sep='\t', header=TRUE, check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
            batch <- batch[sample_ids, args$batch_index]
            names(batch) <- sample_ids
        }else{
            stop('argument --batch is required for: ', args$method)
        }
    }
    # batch removal
    mat <- remove_batch(mat, method=args$method, class_info=class_info, batch_info=batch, ruv_k=args$ruv_k)
    
}  else{
    stop('unknown step: ', args$step)
}   

# write output matrix
message('write expression matrix: ', args$output)
write_matrix(mat, args$output)
