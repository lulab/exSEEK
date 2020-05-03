#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Call domains with significant coverage')
parser$add_argument('-i', '--matrix', type='character', required=TRUE,
    help='input count matrix. Rows are genes. Columns are samples.')
parser$add_argument('-b', '--batch', type='character', required=TRUE,
    help='file containing sample batch information')
parser$add_argument('-o', '--output-file', type='character', required=TRUE,
    help='output file')
parser$add_argument('--batch-index', type='integer', required=FALSE, default=1,
    help='columns index of the batch to be considered')
args <- parser$parse_args()

suppressPackageStartupMessages(library('kBET'))
message('read input batch info: ', args$batch)
batch <-read.table(args$batch, sep='\t', row.names=1, header=TRUE, check.names = FALSE)
message('read input matrix: ', args$matrix)
matrix <- t(read.table(args$matrix, row.names=1, header=TRUE, check.names = FALSE))
batch <- as.data.frame(batch[rownames(matrix),])
batch.estimate <- kBET(matrix, factor(batch[, args$batch_index]), plot=FALSE)
kBET_score <- 1 - mean(batch.estimate$stats$kBET.observed)
message('kBET score: ', kBET_score)
write(kBET_score, args$output_file)
