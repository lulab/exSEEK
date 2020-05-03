#! /usr/bin/env Rscript

suppressPackageStartupMessages(library(argparse))
parser <- ArgumentParser(description='Differential expression')
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
parser$add_argument('-b', '--batch', type='character', required=FALSE,
    help='batch information to remove')
parser$add_argument('--batch-index', type='integer', default=1,
    help='column number of the batch to remove')
parser$add_argument('-m', '--method', type='character', default="deseq2",
    choices=c('deseq2', 'edger_exact', 'edger_glmqlf', 'edger_glmlrt', 'wilcox', 'limma', 'ttest'),
    help='differential expression method to use')
parser$add_argument('--norm-method', type='character', default='TMM',
    choices=c('RLE', 'CPM', 'TMM', 'upperquartile'),
    help='normalization method for count-based methods')
parser$add_argument('--pseudo-count', type='double', default=1.0,
    help='pseudo-count added to log2 transform in ttest')
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

group <- class_info[samples]
mat <- as.matrix(mat[,samples])
#class_info <- as.matrix(class_info)
#colnames(class_info) <- 'label'
mat <- as.matrix(mat)

# read batch information
if(!is.null(args$batch)){
    message('read batch information from: ', args$batch)
    batch <- read.table(args$batch, check.names=FALSE, header=TRUE, as.is=TRUE, row.names=1, sep='\t')
    if((args$batch_index < 1) || (args$batch_index > ncol(batch))){
        stop('Batch index out of bound')
    }
    batch <- batch[samples, args$batch_index]
}else{
    batch <- NULL
}

# set normalization method
if(args$norm_method == 'CPM'){
    norm_method <- 'none'
} else{
    norm_method <- args$norm_method
}
message('perform differential expression using ', args$method)
# Required columns for a differential expression file: baseMean, log2FoldChange, pvalue, padj
if(args$method == 'deseq2'){
    suppressPackageStartupMessages(library(DESeq2))
    if(!is.null(batch)){
        # include batch into regression
        dds <- DESeqDataSetFromMatrix(countData = mat,
                                    colData = as.matrix(data.frame(group=group, batch=batch)),
                                    design = ~group + batch)
    } else {
        dds <- DESeqDataSetFromMatrix(countData = mat,
                                    colData = as.matrix(data.frame(group=group)),
                                    design = ~group)
    }
    dds <- DESeq(dds)
    res <- results(dds, contrast=c('group', 'positive', 'negative'))
    #res <- res[order(res$padj)]
    write.table(as.data.frame(res), args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else if(grepl('^edger_', args$method)) {
    suppressPackageStartupMessages(library(edgeR))
    y <- DGEList(counts=mat, samples=samples, group=group)
    y <- calcNormFactors(y, method=norm_method)
    if(!is.null(batch)){
        # regress out batch information as an additive term
        design <- model.matrix(~group + batch)
    } else {
        design <- model.matrix(~group)
    }
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
        if(!is.null(batch)) message('ignoring batch information for exact text')
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
    matrix_cpm <- cpm(mat, method=norm_method)
    test_func <- function(x){
        wilcox.test(x[group == 'negative'], x[group == 'positive'], alternative='two.sided')$p.value
    }
    pvalues <- apply(matrix_cpm, 1, test_func)
    matrix_logcpm = log2(matrix + args$pseudo_count)
    logFC <- apply(matrix_logcpm[,which(group == 'positive')], 1, mean) -
        apply(matrix_logcpm[,which(group == 'negative')], 1, mean)
    res <- data.frame(log2FoldChange=logFC,
        pvalue=pvalues, 
        padj=p.adjust(pvalues, method='BH'),
        baseMean=apply(matrix_cpm, 1, mean))
    message('Write results to output file: ', args$output_file)
    write.table(res, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else if(args$method == 'limma'){
    suppressPackageStartupMessages(library(limma))
    suppressPackageStartupMessages(library(edgeR))
    y <- DGEList(counts=mat, samples=samples, group=group)
    y <- calcNormFactors(y, method=norm_method)
    if(!is.null(batch)){
        model <- model.matrix(~group + batch)
    } else {
        model <- model.matrix(~group)
    }
    y <- voom(y, model, plot=FALSE)
    fit <- lmFit(y, model)
    fit <- eBayes(fit, robust=TRUE, trend=TRUE)
    #fit2 <- contrasts.ft(fit)
    #fit2 <- eBayes(fit2, robust=TRUE, trend=TRUE)
    #top_table <- topTable(fit2, sort.by='none', n=Inf)
    top_table <- topTable(fit, coef=2, sort.by='none', n=Inf)
    # rename columns
    mapped_names <- colnames(top_table)
    for(i in 1:ncol(top_table)){
        if(colnames(top_table)[i] == 'logFC'){
            mapped_names[i] <- 'log2FoldChange'
        }else if(colnames(top_table)[i] == 'P.Value'){
            mapped_names[i] <- 'pvalue'
        }else if(colnames(top_table)[i] == 'adj.P.Val') {
            mapped_names[i] <- 'padj'
        }else{
            mapped_names[i] <- colnames(top_table)[i]
        }
    }
    colnames(top_table) <- mapped_names

    # write results to file
    message('Write results to output file: ', args$output_file)
    write.table(top_table, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else if(args$method == "ttest") {
    suppressPackageStartupMessages(library(genefilter))
    #mat <- log2(mat + args$pseudo_count)
    res <- rowttests(mat, as.factor(group))
    res$padj <- p.adjust(res$p.value, method='BH')
    res$log2FoldChange <- rowMeans(mat[, group == 'positive']) - rowMeans(mat[, group == 'negative'])
    write.table(res, args$output_file, sep='\t', quote=FALSE, row.names=TRUE)
} else {
    stop('unknown differential expression method: ', args$method)
}
