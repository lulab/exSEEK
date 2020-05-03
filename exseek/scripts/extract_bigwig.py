#! /usr/bin/env python
from __future__ import print_function
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

@command_handler
def abundant_rna_coverage_matrix(args):
    import pandas as pd
    import numpy as np
    from tqdm import tqdm
    import h5py
    #from bx.bbi.bigwig_file import BigWigFile
    from pykent import BigWigFile

    logger.info('read count matrix: ' + args.matrix)
    count_matrix = pd.read_table(args.matrix, sep='\t', index_col=0)
    sample_ids = count_matrix.columns.values
    count_matrix += 1
    read_depth = count_matrix.sum(axis=0)
    cpm = 1e6*count_matrix.div(read_depth, axis=1)
    cpm_mean = np.mean(np.log2(cpm), axis=1)
    # remove miRNA, piRNA, tRNA and genomic RNA
    feature_info = cpm_mean.index.to_series().str.split('|', expand=True)
    feature_info.columns = ['gene_id', 'gene_type', 'gene_name', 'domain_id', 'transcript_id', 'start', 'end']
    feature_info['start'] = feature_info['start'].astype('int')
    feature_info['end'] = feature_info['end'].astype('int')
    cpm_mean = cpm_mean.loc[~feature_info['gene_type'].isin(['piRNA', 'miRNA', 'tRNA', 'genomic'])].copy()
    # sort RNAs by mean log CPM
    order = np.argsort(-cpm_mean)
    # get top n abundant RNAs
    cpm_mean = cpm_mean.iloc[order].head(args.n)
    feature_info = feature_info.loc[cpm_mean.index.values]
    
    # initialize matrix
    matrices = {}
    for feature_name, feature in feature_info.iterrows():
        matrices[feature_name] = np.zeros((len(sample_ids), feature['end']), dtype=np.float32)
    # read bigwig files
    logger.info('read coverage from bigwig files')
    for i_sample, sample_id in tqdm(enumerate(sample_ids), total=len(sample_ids), unit='sample'):
        bigwigs  = {}
        for feature_name, feature in feature_info.iterrows():
            if feature['gene_type'] not in bigwigs:
                #bigwigs[feature['gene_type']] = BigWigFile(open(args.bigwig_pattern.format(sample_id=sample_id, gene_type=feature['gene_type']), 'rb'))
                bigwigs[feature['gene_type']] = BigWigFile(args.bigwig_pattern.format(sample_id=sample_id, gene_type=feature['gene_type']))
            #values = bigwigs[feature['gene_type']].get_as_array(
            #    str(feature['transcript_id']), feature['start'], feature['end'])
            values = bigwigs[feature['gene_type']].query(str(feature['transcript_id']))
            if values is not None:
                matrices[feature_name][i_sample, ~np.isnan(values)] = 1e6*values[~np.isnan(values)]/read_depth[sample_id]
        del bigwigs
    logger.info('write matrices to file: ' + args.output_file)
    with h5py.File(args.output_file, 'w') as fout:
        for feature_name, matrix in tqdm(matrices.items(), unit='matrix'):
            fout.create_dataset(feature_name, data=matrix, compression='gzip')
    fout.close()
    
if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Count reads in BAM files')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('abundant_rna_coverage_matrix',
        help='Extract coverage from abundant long RNAs and build a matrix')
    parser.add_argument('--matrix', type=str, required=True,
        help='count matrix file')
    parser.add_argument('--bigwig-pattern', type=str, required=True,
        help='string format template for bigwig names.'
        'e.g. output/dataset/tbigwig/{sample_id}.{gene_type}.bigWig')
    parser.add_argument('-n', type=int, default=10, 
        help='number of abundant genes to extract')
    parser.add_argument('--output-file', '-o', type=str, required=True,
        help='output matrices in HDF5 format')

    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('extract_bigwig.' + args.command)

    command_handlers.get(args.command)(args)