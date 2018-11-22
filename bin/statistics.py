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
def read_length_hist(args):
    import pysam
    import numpy as np
    from ioutils import open_file_or_stdout

    logger.info('read input BAM/SAM file: ' + args.input_file)
    sam = pysam.AlignmentFile(args.input_file, "rb")
    counts_ref = np.zeros(args.max_length, dtype=np.int64)
    counts_query = np.zeros(args.max_length, dtype=np.int64)
    max_length = args.max_length
    for read in sam:
        counts_query[read.query_length] += 1
        counts_ref[min(read.reference_length, max_length - 1)] += 1
    
    logger.info('create output file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as f:
        f.write('length\tquery\treference\n')
        for i in range(args.max_length):
            f.write('{}\t{}\t{}\n'.format(i, counts_query[i], counts_ref[i]))

@command_handler
def read_duplicate_hist(args):
    import pysam
    import numpy as np
    from ioutils import open_file_or_stdout

    bin_size = args.bin_size
    max_length = args.max_length
    n = max_length//bin_size
    bounds = np.arange(0, (n + 1)*bin_size, bin_size)

    logger.info('read chrom sizes: ' + args.chrom_sizes_file)
    chrom_sizes = {}
    with open(args.chrom_sizes_file, 'r') as f:
        for line in f:
            c = line.strip().split('\t')
            chrom_sizes[c[0]] = int(c[1])
    logger.info('read input BAM/SAM file: ' + args.input_file)
    sam = pysam.AlignmentFile(args.input_file, "rb")
    dup_counts = np.zeros(n + 1, dtype=np.int64)
    tot_counts = np.zeros(n + 1, dtype=np.int64)
    max_length = args.max_length
    for read in sam:
        index = min(chrom_sizes[read.reference_name]//bin_size, n)
        if read.is_duplicate:
            dup_counts[index] += 1
        tot_counts[index] += 1
    
    logger.info('create output file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as f:
        f.write('bin\tduplicates\ttotal\n')
        for i in range(n + 1):
            f.write('{}\t{}\t{}\n'.format(bounds[i], dup_counts[i], tot_counts[i]))

@command_handler
def fragment_length_hist(args):
    import pysam
    import numpy as np
    from ioutils import open_file_or_stdout

    logger.info('read input BAM/SAM file: ' + args.input_file)
    sam = pysam.AlignmentFile(args.input_file, "rb")
    max_length = args.max_length
    counts = np.zeros(max_length + 1, dtype=np.int64)
    read1 = None
    for read in sam:
        if (not read.is_paired) or (not read.is_proper_pair):
            continue
        if read.is_read1:
            read1 = read
        elif read.is_read2:
            if read.query_name != read1.query_name:
                continue
            length = read.reference_end - read1.reference_start
            counts[min(length, max_length)] += 1
    
    with open_file_or_stdout(args.output_file) as f:
        f.write('fragment_length\tcounts\n')
        for i in range(max_length + 1):
            f.write('{}\t{}\n'.format(i, counts[i]))
        

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Statistics of exRNA datasets')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('read_length_hist', 
        help='calculate a histogram for read length distribution')
    parser.add_argument('--input-file', '-i', type=str, required=True, help='input BAM/SAM file')
    parser.add_argument('--output-file', '-o', type=str, default='-', help='output histogram file')
    parser.add_argument('--max-length', '-l', type=int, default=10000)

    parser = subparsers.add_parser('read_duplicate_hist', 
        help='calculate a histogram for duplicate reads fraction vs read length')
    parser.add_argument('--input-file', '-i', type=str, required=True, help='input BAM/SAM file')
    parser.add_argument('--output-file', '-o', type=str, default='-', help='output histogram file')
    parser.add_argument('--chrom-sizes-file', type=str, required=True, 
        help='file containing chromosome sizes')
    parser.add_argument('--bin-size', type=int, default=10, 
        help='bin size for read length')
    parser.add_argument('--max-length', type=int, default=10000,
        help='upper bound for bins')

    parser = subparsers.add_parser('fragment_length_hist',
        help='calculate a histogram for fragment length in a paired-end BAM/SAM file')
    parser.add_argument('--input-file', '-i', type=str, required=True, help='input BAM/SAM file')
    parser.add_argument('--output-file', '-o', type=str, default='-', help='output histogram file')
    parser.add_argument('--max-length', type=int, default=1000,
        help='upper bound for fragment lengths')
    
    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('statistics.' + args.command)

    command_handlers.get(args.command)(args)