#!/usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()
parser$add_argument("-s", "--step", required=TRUE, default='imputation', help="which step to run")
parser$add_argument("-i", "--input", required=TRUE, help="input expression matrix file")
parser$add_argument("-c", "--class", required=TRUE, help="input class info file")
parser$add_argument("-b", "--batch", required=TRUE, help="input batch info file")

parser$add_argument("--filterout", required=TRUE, help="output filter path")
parser$add_argument("--imputeout", required=TRUE, help="output imputation path")
parser$add_argument("--normalizeout", required=TRUE, help="output normalization file")
parser$add_argument("--batchremoveout", required=TRUE, help="output batchremoved file")

parser$add_argument("--filtercount", type="integer", default=5,
    help="filter by counts of a gene [default = %(default)s]",
    metavar="NUMBER")
parser$add_argument("--filtersample", type="integer", default=10,
    help="filter by counts of sample above certain counts of a gene [default = %(default)s]",
    metavar="NUMBER")
parser$add_argument( "--imputemethod", type="character", default="scimpute_count",
                    metavar="STRING",
                    help="the imputation algorithm to use [default = %(default)s]")
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

parser$add_argument( "--normmethod", type="character", default="SCNorm",
                    metavar="STRING",
                    help="the normalization algorithm to use [default = %(default)s]")
parser$add_argument( "--normtopk", type="integer", default=20,
                    metavar="NUMBER",
                    help="top K feature as scale factor [default = %(default)s]")
parser$add_argument( "--cvthreshold", type="double", default=0.5,
                    metavar="NUMBER",
                    help="coefficient variance threshold of reference gene, filter ref gene with CV bigger than [default = %(default)s]")
parser$add_argument( "--removetype", type="character", default="miRNA,piRNA",
                    metavar="STRING",
                    help="remove some time of RNA for normalization scale factor calculation [default = %(default)s]")
parser$add_argument( "--refergenefile", type="character", #default="miRNA,piRNA",
                    metavar="STRING",
                    help="reference gene file path [default = %(default)s]")
#they are feature name for full length feature, most are miRNA, for domain feature, they have the same feature name

parser$add_argument("--batchmethod", type="character", default="RUV",
                    metavar="STRING",
                    help="the batch removal algorithm to use [default = %(default)s]")
parser$add_argument("--batchindex", type="integer", default=1,
                    metavar="INT",
                    help="batch index to select which batch to use [default = %(default)s]")

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
scimpute_count <- function(mat,tmp_path=".",impute_path=".", K = 5, N = 3) {
    suppressMessages(library("scImpute"))
    print('start imputation using scImpute')
    mat_correct <- names(mat)
    names(mat) <-  paste('C_',seq_len(length(names(mat))))
    write.csv(mat, paste(tmp_path,"tmpsave.csv",sep=""), sep=',')
    scimpute(count_path = paste(tmp_path,"tmpsave.csv",sep=""), infile = "csv", outfile = "txt", out_dir = impute_path , Kcluster = K, ncores = N)
    mat <- read.table(paste(impute_path,"scimpute_count.txt",sep=""),sep=' ',header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
    names(mat) <-mat_correct
    write.table(  mat , file=paste(args$imputeout,'filter.scimpute_count.',splitname,sep='') ,sep='\t' )
}

viper_count <- function(mat,tmp_path=".",impute_path=".",num = 5000, percentage.cutoff = 0.1, alpha= 0.5) {
    suppressWarnings(library(VIPER))
    mat_correct <- names(mat)
    names(mat) <-  paste('C_',seq_len(length(names(mat))))
    print (num, percentage.cutoff, alpha)
    VIPER(mat, num = num, percentage.cutoff = percentage.cutoff, minbool = FALSE, alpha = alpha, report = TRUE, outdir = impute_path)
    mat <- read.table(paste(impute_path,'imputed_counts.csv',sep=''),sep=' ',header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
    names(mat) <-mat_correct
    write.table(mat, file=paste(args$imputeout,'filter.viper_count.',splitname,sep='') ,sep='\t' )
}


################################################################################
###############################normalization####################################
################################################################################
library(magrittr)
normalize_check_arg <- function(norm_methods, top_n, rm_gene_type, refer_gene_id) {
    if ('CPM_top' %in% norm_methods && is.null(top_n))
    stop('top_n must be specified for CPM_top method')
    
    if ('CPM_rm' %in% norm_methods && is.null(rm_gene_type))
    stop('rm_gene_type must be specified for CPM_rm method')
    
    if ('CPM_refer' %in% norm_methods && is.null(refer_gene_id))
    stop('refer_gene_id must be specified for CPM_refer method')
    
}

#' @title martix normalization
#' @examples
#' \donotrun{
#'     norm_mat(
#'         '/path/to/matrix'
#'     )
#' }
normalize <- function(
mat,
norm_methods = c( 'SCnorm', 'TMM', 'RLE', 'CPM', 'CPM_top', 'CPM_rm', 'CPM_refer', 'null'),
top_n = 20, rm_gene_type = 'miRNA,tRNA',
sample_class_path = NULL, PCA_label_by = NULL, PAC_color_by = NULL,
refer_gene_id_path='data/matrix_processing/refer_gene_id.txt',
tmp_path='.',impute_path="./imputation/", K = 5, N = 3,
output_dir = '.',output_file = 'norm',cv_threshold=1,imputemethod = 'scimpute_count'
) {
    rm_gene_type <- unlist(strsplit(rm_gene_type,',', fixed = TRUE))
    normalize_check_arg(norm_methods, top_n, rm_gene_type, read.table(refer_gene_id_path,header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)[,1])
    if ('SCnorm' %in% norm_methods)    Norm_SCnorm <- norm_SCnorm(mat)
    if ('TMM' %in% norm_methods)       Norm_TMM <- norm_tmm(mat)
    if ('RLE' %in% norm_methods)       Norm_RLE <- norm_rle(mat)
    if ('CPM' %in% norm_methods)       Norm_CPM <- norm_cpm_total(mat)
    if ('CPM_top' %in% norm_methods)   Norm_CPM_top <- norm_cpm_top(mat, top_n)
    if ('CPM_rm' %in% norm_methods)    Norm_CPM_rm <- norm_cpm_rm(mat, rm_gene_type)
    if ('CPM_refer' %in% norm_methods) Norm_CPM_refer <- norm_cpm_refer(mat, refer_gene_id_path)
    if ('null' %in% norm_methods)      Norm_null <- mat

    
    tmp_names <- ls(pattern = '^Norm_')
    imputename <-paste(imputemethod,'.',sep='')
    for(tmp_name in tmp_names) {
        write.table(get(tmp_name), paste(output_dir,'filter.',imputename, tmp_name,'.',splitname,sep=''),sep='\t')
    }
#    if (imputemethod=='null'){
#        write.table(mat, paste(output_dir,'filter.',imputename,'Norm_null.',splitname,sep=''),sep='\t')
#    }
    
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
    sce <- suppressMessages(SCnorm::SCnorm(mat, Conditions, ...));
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
#' @details `norm_tmm()` performs TMM normalization
#'
#' @examples
#' norm_tmm(sim_mat)
#'
#' @export
norm_tmm <- function(mat) {
    print('start normalization using TMM')
    mat %>% as_SingleCellExperiment() %>%
    {suppressWarnings(scater::normaliseExprs(., "TMM"))} %>%
    scater::normalise() %>% SingleCellExperiment::normcounts()
}


#' @rdname  norm_scater
#'
#' @details `norm_rle()` performs RLE normalization
#'
#' @examples
#' norm_rle(sim_mat)
#'
#' @export
norm_rle <- function(mat) {
    print('start normalization using RLE')
    mat %>% as_SingleCellExperiment() %>%
    {suppressWarnings(scater::normaliseExprs(., "RLE"))} %>%
    scater::normalise() %>% SingleCellExperiment::normcounts()
}


# norm_cpm ------------------

#' @title CPM normalization by some genes
#'
#' @param mat integer matrix. counts
#' @param row integer or logical. Use which rows (genes) as normalization factor
norm_cpm_impl <- function(mat, row) {
    t(t(mat*1e6) / colSums(mat[row, , drop = F], na.rm = T))
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


norm_cpm_refer <- function(mat, refer_gene_id_path='data/matrix_processing/refer_gene_id.txt',cv_threshold=0.5) {
    print(paste('start normalization by reference gene ','thresholding reference gene by Coefficient of Variance:',args$cvthreshold,sep=' '))
    refer_gene_id <- read.table(refer_gene_id_path,header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)[,1]
    row_refer <- mat %>% rownames() %>% stringr::str_extract('[^|]+') %>%
    {. %in% refer_gene_id}
    if (!any(row_refer))
    stop('can\'t find any reference genes in the matrix for CPM normalization')
    keeped_ref <- mat[row_refer, , drop = F] %>% apply(1, cv_fun) < cv_threshold
    norm_cpm_impl(mat, row_refer[keeped_ref])
}

################################################################################
#################################batch removal##################################
################################################################################


batch <- function(
classinfo_path,
batchinfo_path,
input_path,
norm_methods = 'SCnorm',
imputemethod = 'scimpute_count',
batch_methods = c('ruv','combat'),
output_path = './',batchindex
){
    suppressMessages(library(EDASeq))
    suppressMessages(library(RUVSeq))
    suppressMessages(library(sva))
    suppressMessages(library(scRNA.seq.funcs))
    
    if ('SCnorm' %in% norm_methods)    normname <- 'Norm_SCnorm'
    if ('TMM' %in% norm_methods)       normname <- 'Norm_TMM'
    if ('RLE' %in% norm_methods)       normname <- 'Norm_RLE'
    if ('CPM' %in% norm_methods)       normname <- 'Norm_CPM'
    if ('CPM_top' %in% norm_methods)   normname <- 'Norm_CPM_top'
    if ('CPM_rm' %in% norm_methods)    normname <- 'Norm_CPM_rm'
    if ('CPM_refer' %in% norm_methods) normname <- 'Norm_CPM_refer'
    if ('null' %in% norm_methods) normname <- 'Norm_null'
    #mat <- get(paste('mat_',tolower(norm_methods),sep=''))
    
    
    imputename <-paste(imputemethod,'.',sep='')
    normname <-paste(normname,'.',sep='')

    mat <- read.table(paste(input_path,'filter.',imputename, normname,splitname,sep=''),sep='\t',header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
    print (dim(mat))
    if ('RUV' %in% batch_methods)       Batch_RUV <- ruv(mat,classinfo_path,output_path,2,10)
    if ('Combat' %in% batch_methods)    Batch_Combat <- combat(mat,batchinfo_path,output_path,batchindex)
    tmp_names <- ls(pattern = '^Batch_')
    for(tmp_name in tmp_names) {
        if (tmp_name=='Batch_Combat'){
        write.table(get(tmp_name), paste(output_path,'filter.',imputename, normname, tmp_name,'_',toString(batchindex),'.',splitname,sep=''),sep='\t',quote=FALSE, row.names=TRUE, col.names=TRUE)
        } else {
        write.table(get(tmp_name), paste(output_path,'filter.',imputename, normname, tmp_name,'.',splitname,sep=''),sep='\t',quote=FALSE, row.names=TRUE, col.names=TRUE)
        }
    }
    if (batch_methods=='null'){
        write.table(mat, paste(output_path,'filter.',imputename, normname,'Batch_null.',splitname,sep=''),sep='\t',quote=FALSE, row.names=TRUE, col.names=TRUE)
    }
}
ruv <- function(
    mat,
    classinfo_path,
    output_path,
    label_column = 2,
    k = 10
){
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
    #write.table(ruv$normalizedCounts, file=paste(output_path,'mx_ruv_batchremoval.txt',sep=''), sep='\t', quote=FALSE, row.names=TRUE, col.names=TRUE)
    exp(ruv$normalizedCounts)
}


combat <- function(
    mat,
    batchinfo_path,
    output_path,
    batch_column = batchindex
){
    print('start batch removal using combat')
    batch_info <-read.table(batchinfo_path,sep='\t',row.names=1,header=T,check.names = FALSE)
    batchname <-toString(names(batch_info)[batch_column])
    batch_info=batch_info[names(mat),]
    mod <- model.matrix(~ 1, data = batch_info)
    if (!dim(mat)[2]==dim(batch_info)[1])
    stop('sample numbers in batch info and expression matrix should be same')
    combat <- ComBat(
        dat = log(mat+0.001),
        batch = factor(batch_info[,batch_column]),
        mod = mod,
        par.prior = TRUE,
        prior.plots = FALSE
    )
    
    mat <- exp(combat)
    #write.csv(mat, file=paste(output_path,'mx_combat_batchremoval_',batchname,'.txt',sep=''))
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


# filter
if (args$step =='filter'){
mat_raw <- read.table(args$input,sep='\t',header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
mat_filter <-filter_low(mat_raw,args$filtercount, args$filtersample)
write.table(mat_filter, paste(args$filterout, 'filter','.',splitname,sep=''),sep='\t')
} else if(args$step =='imputation'){
# imputation
library(readr)
mat_filter <-read.table(paste(args$filterout,'filter.',splitname,sep=''),sep='\t',header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
if (args$imputemethod=='scimpute_count'){
scimpute_count(mat_filter,impute_path= args$imputeout, K = args$imputecluster, N = args$processors)
} else if(args$imputemethod=='viper_count'){
viper_count(mat_filter, impute_path= args$imputeout,num = args$imputevipernum,percentage.cutoff = args$imputecutoff, alpha= args$imputealpha) 
} else if(args$imputemethod=='null'){
write.table( mat_filter, paste(args$imputeout,'filter.null.',splitname,sep='') ,sep='\t')
}
}  else if(args$step =='normalization'){
# normalization
imputemethod <-args$imputemethod
imputename <-paste(imputemethod,'.',sep='')
mat_impute <- read.table(paste(args$imputeout,'filter.',imputename,splitname,sep=""),header=TRUE,  check.names=FALSE, row.names=1, stringsAsFactors=FALSE)
normalize(mat_impute, norm_methods =args$normmethod,top_n = args$normtopk,cv_threshold = args$cvthreshold, imputemethod=args$imputemethod,
rm_gene_type = args$removetype,refer_gene_id_path = args$refergenefile, output_dir = args$normalizeout, K = args$imputecluster, N = args$processors)
}  else if(args$step =='batch_removal'){
# batch removal
batch(classinfo_path = args$class, batchinfo_path = args$batch,
output_path=args$batchremoveout,input_path = args$normalizeout,imputemethod=args$imputemethod,
norm_methods = args$normmethod,batch_methods = args$batchmethod,args$batchindex)
}  else{
    print ('unexpected step')
}   
