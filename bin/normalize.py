#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

def normalize(args):
    from ioutils import open_file_or_stdin, open_file_or_stdout
    import pandas as pd

    with open_file_or_stdin(args.input_file) as f:
        matrix = pd.read_table(f, sep='\t', index_col=0)
    if args.method == 'cpm':
        matrix = 1e6*matrix.astype('float')/matrix.sum(axis=0)
    with open_file_or_stdout(args.output_file) as f:
        matrix.to_csv(f, sep='\t', header=True, index=True, na_rep='NA')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Normalization module')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='input feature matrix (rows are samples and columns are features')
    parser.add_argument('--method', '-m', type=str, default='cpm',
        choices=('cpm',), help='normalization method')
    parser.add_argument('--transpose', '-t', action='store_true', 
        help='transpose the matrix before normalization')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='normalized matrix file')
    
    args = parser.parse_args()
    logger = logging.getLogger('normalize')
    normalize(args)
