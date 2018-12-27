#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Call domains with significant coverage')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-c', '--classes', type='character', required=TRUE,
    help='input sample class information. Column 1: sample_id, Column 2: class')
parser$add_argument('-p', '--positive-class', type='character', required=TRUE,
    help='comma-separated class names to use as positive class')
parser$add_argument('-n', '--negative-class', type='character', required=TRUE,
    help='comma-separated class names to use as negative class')
parser$add_argument('-m', '--method', type='character', default="deseq2",
    choices=c('deseq2'), help='differential expression method to use')
parser$add_argument('-o', '--output-file', type='character', required=TRUE,
    help='output file')
args <- parser$parse_args()

message('read count matrix: ', args$matrix)
mat <- read.table(args$matrix, header = TRUE, row.names=1, check.names=FALSE, sep='\t')
message('read class information: ', args$classes)
class_info <- read.table(args$classes, row.names=1, check.names=FALSE, header = TRUE, sep='\t', as.is=TRUE)
class_info <- class_info[colnames(mat),]
names(class_info) <- colnames(mat)
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

if(args$method == 'deseq2'){
    suppressPackageStartupMessages(library(DESeq2))
    dds <- DESeqDataSetFromMatrix(countData = mat,
                                colData = class_info,
                                design = ~label
                                    )
    dds <- DESeq(dds)
    res <- results(dds, contrast=c('label', 'positive', 'negative'))
    #res <- res[order(res$padj)]
    write.table(as.data.frame(res), args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else {
    stop('unknown differential expression method: ', args$method)
}
