#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Call domains with significant coverage')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-c', '--classes', type='character', required=TRUE,
    help='input sample class information. Column 1: sample_id, Column 2: class')
parser$add_argument('-s', '--samples', type='character',
    help='input file containing sample ids for differential expression')
parser$add_argument('-p', '--positive-class', type='character', required=TRUE,
    help='comma-separated class names to use as positive class')
parser$add_argument('-n', '--negative-class', type='character', required=TRUE,
    help='comma-separated class names to use as negative class')
parser$add_argument('-m', '--method', type='character', default="deseq2",
    choices=c('deseq2', 'edger_exact', 'edger_glmqlf', 'edger_glmlrt', 'wilcox'), 
    help='differential expression method to use')
parser$add_argument('-o', '--output-file', type='character', required=TRUE,
    help='output file')
args <- parser$parse_args()

message('read count matrix: ', args$matrix)
mat <- read.table(args$matrix, header = TRUE, row.names=1, check.names=FALSE, sep='\t')
message('read class information: ', args$classes)
class_info <- read.table(args$classes, row.names=1, check.names=FALSE, header = TRUE, sep='\t', as.is=TRUE)
class_info <- class_info[colnames(mat),]
names(class_info) <- colnames(mat)
if(!is.null(args$samples)){
    message('read sample ids: ', args$samples)
    samples <- read.table(args$samples, check.names=FALSE, header=FALSE, as.is=TRUE)[,1]
    mat <- mat[, samples]
    class_info <- class_info[samples]
}
# get positive and negative class
for(cls in strsplit(args$positive_class, ',')[[1]]){
    class_info[class_info == cls] <- 'positive'
}
positive_samples <- names(class_info)[class_info == 'positive']
if(length(positive_samples) == 0){
    stop('No positive samples found')
}
message('Number of positive samples: ', length(positive_samples))
negative_samples <- NULL
for(cls in strsplit(args$negative_class, ',')[[1]]){
    class_info[class_info == cls] <- 'negative'
}
negative_samples <- names(class_info)[class_info == 'negative']
if(length(negative_samples) == 0){
    stop('No negative samples found')
}
message('Number of negative samples: ', length(negative_samples))

samples <- c(positive_samples, negative_samples)

class_info <- class_info[samples]
mat <- mat[,samples]
class_info <- as.matrix(class_info)
colnames(class_info) <- 'label'
mat <- as.matrix(mat)

# Required columns for a differential expression file: baseMean, log2FoldChange, pvalue, padj
if(args$method == 'deseq2'){
    suppressPackageStartupMessages(library(DESeq2))
    dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData = class_info,
                                design = ~label)
    dds <- DESeq(dds)
    res <- results(dds, contrast=c('label', 'positive', 'negative'))
    #res <- res[order(res$padj)]
    write.table(as.data.frame(res), args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else if(grepl('^edger_', args$method)) {
    suppressPackageStartupMessages(library(edgeR))
    group <- class_info[,1]
    y <- DGEList(counts=mat, samples=samples, group=group)
    y <- calcNormFactors(y)
    design <- model.matrix(~group)
    y <- estimateDisp(y, design)
    if(args$method == 'edger_glmqlf'){
        fit <- glmQLFit(y, design)
        test <- glmQLFTest(fit, coef=2)
        res <- topTags(test, n=nrow(mat), sort.by='none')
    } else if(args$method == 'edger_glmlrt'){
        fit <- glmFit(y, design)
        test <- glmLRT(fit, coef=2)
        res <- topTags(test, n=nrow(mat), sort.by='none')
    } else if(args$method == 'edger_exact'){
        test <- exactTest(y)
        res <- topTags(test, n=nrow(mat), sort.by='none')
    }
    res <- cbind(res$table, baseMean=2^(res$table$logCPM))
    # rename columns
    mapped_names <- colnames(res)
    for(i in 1:ncol(res)){
        if(colnames(res)[i] == 'logFC'){
            mapped_names[i] <- 'log2FoldChange'
        }else if(colnames(res)[i] == 'PValue'){
            mapped_names[i] <- 'pvalue'
        }else if(colnames(res)[i] == 'FDR') {
            mapped_names[i] <- 'padj'
        }else{
            mapped_names[i] <- colnames(res)[i]
        }
    }
    colnames(res) <- mapped_names

    # write results to file
    message('Write results to output file: ', args$output_file)
    write.table(res, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else if(args$method == 'wilcox') {
    suppressPackageStartupMessages(library(edgeR))
    # normalize
    matrix_cpm <- cpm(mat)
    group <- class_info[,1]
    test_func <- function(x){
        wilcox.test(x[group == 'negative'], x[group == 'positive'], alternative='two.sided')$p.value
    }
    pvalues <- apply(matrix_cpm, 1, test_func)
    logFC <- apply(log2(matrix_cpm[,which(group == 'positive')]), 1, mean) -
        apply(log2(matrix_cpm[,which(group == 'negative')]), 1, mean)
    res <- data.frame(log2FoldChange=logFC,
        pvalue=pvalues, 
        padj=p.adjust(pvalues, method='BH'),
        baseMean=apply(matrix_cpm, 1, mean))
    message('Write results to output file: ', args$output_file)
    write.table(res, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else {
    stop('unknown differential expression method: ', args$method)
}
