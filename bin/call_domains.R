#! /usr/bin/env Rscript
library(countreg)

args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 2) {
    message("Usage: call_domains.R input_file fdr")
    quit(save='no', status=1)
}
bed_file <- args[1]
fdr_value <- as.numeric(args[2])

message('read coverage file:', bed_file)
bed <- read.table(bed_file, sep='\t', header=FALSE)
message('numbe of bins:', nrow(bed))
# remove bins with 1 read
bed <- bed[bed[, 5] > 1,]
x <- as.integer(bed[, 5])

ztnb.test <- function(x, mu, theta) {
    # build a p-value table
    max_x = qztnbinom(1 - 1e-10, mu=mu, theta=theta)
    pvalue_table <- dztnbinom(1:max_x, mu=mu, theta=theta)
    pvalue_table <- rev(cumsum(rev(pvalue_table)))
    # calculate pvalues
    pvalues <- rep(1e-10, length(x))
    pvalues[x <= max_x] <- pvalue_table[x[x <= max_x]]
    qvalues <- p.adjust(pvalues, method='BH')
    return(qvalues)
}

bg = x[x < quantile(x, 0.99)]
message('fit ZTNB distribution')
fit <- zerotrunc(bg ~ 1, dist='negbin')
qvalues <- ztnb.test(x, mu=exp(fit$coefficients), theta=fit$theta)
bg <- x[qvalues >= fdr_value]
message('fit ZTNB distribution')
fit <- zerotrunc(bg ~ 1, dist='negbin')
qvalues <- ztnb.test(x, mu=exp(fit$coefficients), theta=fit$theta)
bed_sig <- bed[qvalues < fdr_value,]
bed_sig <- cbind(bed_sig, qvalues[qvalues < fdr_value])
message('number of significant bins:', nrow(bed_sig))
message('write significant domains')
write.table(bed_sig, row.names=FALSE, col.names=FALSE, sep='\t', quote=FALSE)
