#! /usr/bin/env Rscript
suppressPackageStartupMessages(library("argparse"))
parser <- ArgumentParser(description='Call domains with significant coverage')
parser$add_argument('-p', '--pvalue', type='double', default=0.05,
    help='adjusted p-value threshold for defining significant domains', metavar='NUMBER')
parser$add_argument('-i', '--input-file', type='character', default='-',
    help='input BED file with covearge in the 5th column')
parser$add_argument('-o', '--output-file', type='character', default='-',
    help='output BED file with adjusted p-value in the last column')
parser$add_argument('-b', '--background', type='double', default=0.99,
    help='set coverage quantile below this quantile as background regions')
args <- parser$parse_args()

library(countreg)

ztnb.test <- function(x, mu, theta, max_x=2000) {
    # build a p-value table
    pvalue_table <- dztnbinom(1:max_x, mu=mu, theta=theta)
    pvalue_table <- rev(cumsum(rev(pvalue_table)))
    # calculate pvalues
    pvalues <- rep(pvalue_table[max_x], length(x))
    pvalues[x <= max_x] <- pvalue_table[x[x <= max_x]]
    qvalues <- p.adjust(pvalues, method='BH')
    return(qvalues)
}
message('read coverage file:', args$input_file)
if(args$input_file == '-'){
    bed <- read.table(file('stdin'), sep='\t', header=FALSE)
} else {
    bed <- read.table(args$input_file, sep='\t', header=FALSE)
}
message('number of bins:', nrow(bed))
# remove bins with 1 read
#bed <- bed[bed[, 5] > 1,]
x <- as.integer(bed[, 5])
bg = x[x < quantile(x, args$background)]
if(length(bg) > 0){
    message('fit ZTNB distribution')
    fit <- zerotrunc(bg ~ 1, dist='negbin')
    qvalues <- ztnb.test(x, mu=exp(fit$coefficients), theta=fit$theta)
    bg <- x[qvalues >= args$pvalue]
    message('fit ZTNB distribution')
    fit <- zerotrunc(bg ~ 1, dist='negbin')
    qvalues <- ztnb.test(x, mu=exp(fit$coefficients), theta=fit$theta)
    bed_sig <- bed[qvalues < args$pvalue,]
    bed_sig <- cbind(bed_sig, qvalues[qvalues < args$pvalue])
} else {
    # empty file
    bed_sig <- data.frame()
}
message('number of significant bins:', nrow(bed_sig))
message('write significant domains')
if(args$output_file == '-'){
    write.table(bed_sig, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
} else {
    write.table(bed_sig, args$output_file, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
}
