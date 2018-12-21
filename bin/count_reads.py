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

@command_handler
def count_circrna(args):
    import HTSeq
    import numpy as np
    import pandas as pd
    from collections import OrderedDict
    from ioutils import open_file_or_stdout

    logger.info('read input BAM/SAM file: ' + args.input_file)
    if args.input_file.endswith('.sam'):
        sam = HTSeq.SAM_Reader(args.input_file)
    elif args.input_file.endswith('.bam'):
        sam = HTSeq.BAM_Reader(args.input_file)
    else:
        raise ValueError('unsupported file extension')
    
    # extract junction positions from SAM header
    logger.info('extract junction positions')
    junction_positions = OrderedDict()
    for sq in sam.get_header_dict()['SQ']:
        junction_positions[sq['SN']] = sq['LN']//2
    # initialize counts
    gene_ids = list(junction_positions.keys())
    counts = pd.Series(np.zeros(len(gene_ids), dtype='int'), index=gene_ids)
    # count reads
    min_mapping_quality = args.min_mapping_quality
    if args.paired_end:
        logger.info('count paired-end fragments')
        for bundle in HTSeq.pair_SAM_alignments(sam, bundle=True):
            # ignore multi-mapped pairs
            if len(bundle) != 1:
                continue
            read1, read2 = bundle[0]
            # ignore singletons
            if read1 is None or read2 is None:
                continue
            # ignore unmapped reads
            if not read1.aligned and read2.aligned:
                continue
            # ignore pairs with mapping quality below threshold
            if (read1.aQual < min_mapping_quality) or (read2.aQual < min_mapping_quality):
                continue
            # ignore pairs on different chromosomes
            if read1.iv.chrom != read2.iv.chrom:
                continue
            pos = junction_positions[read1.iv.chrom]
            if read1.iv.start < pos <= read2.iv.end:
                counts[read1.iv.chrom] += 1
    else:
        logger.info('count single-end reads')
        for read in sam:
            # ignore unmapped read
            if not read.aligned:
                continue
            # ignore reads with mapping quality below threshold
            if read.aQual < min_mapping_quality:
                continue
            pos = junction_positions[read.iv.chrom]
            if read.iv.start < pos <= read.iv.end:
                counts[read.iv.chrom] += 1
    # output counts
    logger.info('count fragments: {}'.format(counts.sum()))
    logger.info('write counts to file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as fout:
        counts.to_csv(fout, sep='\t', header=None, index=True, na_rep='NA')
                
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

    parser = subparsers.add_parser('count_circrna', 
        help='count reads/fragments mapped to circRNA junctions')
    parser.add_argument('--input-file', '-i', type=str, required=True, help='input BAM/SAM file')
    parser.add_argument('--paired-end', '-p', action='store_true', help='count reads as paired-end')
    parser.add_argument('--min-mapping-quality', '-q', type=int, default=0,
        help='only count reads with mapping quality greater than this number')
    parser.add_argument('--output-file', '-o', type=str, default='-', 
        help='output tab-deliminated file. Two columns: gene_id, count')
    
    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('count_reads.' + args.command)

    command_handlers.get(args.command)(args)