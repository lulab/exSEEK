#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

@command_handler
def count_transcript(args):
    import pysam
    import numpy as np
    from ioutils import open_file_or_stdout
    from collections import OrderedDict

    logger.info('read input BAM/SAM file: ' + args.input_file)
    sam = pysam.AlignmentFile(args.input_file, "rb")
    counts = OrderedDict()
    min_mapping_quality = args.min_mapping_quality
    strandness = {'no': 0, 'forward': 1, 'reverse': 2}.get(args.strandness, 0)
    for read in sam:
        if read.is_unmapped:
            continue
        if read.mapping_quality < min_mapping_quality:
            continue
        if (strandness == 1) and read.is_reverse:
            continue
        if (strandness == 2) and (not read.is_reverse):
            continue
        if read.reference_name not in counts:
            counts[read.reference_name] = 0
        counts[read.reference_name] += 1
    
    with open_file_or_stdout(args.output_file) as f:
        if sam.header is not None:
            for name in sam.header.references:
                f.write('{}\t{}\n'.format(name, counts.get(name, 0)))
        else:
            for name, count in counts.items():
                f.write('{}\t{}\n'.format(name, count))

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Count reads in BAM files')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('count_transcript', 
        help='count reads in BAM in transcript coordinates')
    parser.add_argument('--input-file', '-i', type=str, required=True, help='input BAM/SAM file')
    parser.add_argument('--min-mapping-quality', '-q', type=int, default=0,
        help='only count reads with mapping quality greater than this number')
    parser.add_argument('--strandness', '-s', type=str, default='no',
        choices=('forward', 'reverse', 'no'),
        help='forward/reverse: only count reads in reverse strand. no: count reads in both strands')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output file')
    
    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('count_reads.' + args.command)

    command_handlers.get(args.command)(args)