#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

@command_handler
def bigwig_to_hdf5(args):
    import h5py
    from pykent import init_library, BigWigFile
    from tqdm import tqdm
    
    logger.info('initialize shared library')
    init_library(args.library)
    logger.info('read input BigWig file: ' + args.input_file)
    bwf = BigWigFile(args.input_file)
    logger.info('open output HDF5 file: ' + args.output_file)
    fout = h5py.File(args.output_file, 'w')
    chrom_sizes = bwf.get_chrom_list()
    for chrom, size in tqdm(chrom_sizes, unit='chrom'):
        values = bwf.values_on_chrom(chrom, fillna=args.fillna)
        fout.create_dataset(chrom, data=values, compression='gzip')
    fout.close()
    
if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Count reads in BAM files')
    main_parser.add_argument('--library', '-l', type=str,
        help='path to libjkweb.so')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('bigwig_to_hdf5')
    parser.add_argument('--input-file', '-i', type=str, required=True,
        help='input BigWig file')
    parser.add_argument('--output-file', '-o', type=str, required=True,
        help='output HDF5 file')
    parser.add_argument('--fillna', type=float)

    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('pyucsctools.' + args.command)

    command_handlers.get(args.command)(args)
