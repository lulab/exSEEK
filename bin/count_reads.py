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
            for sq in sam.header['SQ']:
                name = sq['SN']
                f.write('{}\t{}\n'.format(name, counts.get(name, 0)))
        else:
            for name, count in counts.items():
                f.write('{}\t{}\n'.format(name, count))

@command_handler
def count_circrna(args):
    import HTSeq
    import numpy as np
    import pandas as pd
    from collections import OrderedDict, defaultdict
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
    strandness = args.strandness
    if args.paired_end:
        logger.info('count paired-end fragments')
        stats = defaultdict(int)
        for bundle in HTSeq.pair_SAM_alignments(sam, bundle=True):
            stats['total_pairs'] += 1
            # ignore multi-mapped pairs
            if len(bundle) != 1:
                stats['multi_mapping'] += 1
                continue
            read1, read2 = bundle[0]
            # ignore singletons
            if (read1 is None) or (read2 is None):
                stats['singleton'] += 1
                continue
            # ignore unmapped reads
            if not (read1.aligned and read2.aligned):
                stats['unmapped'] += 1
                continue
            # ignore pairs with mapping quality below threshold
            if (read1.aQual < min_mapping_quality) or (read2.aQual < min_mapping_quality):
                stats['low_mapping_quality'] += 1
                continue
            if (strandness == 'forward') and (not ((read1.iv.strand == '+') and (read2.iv.strand == '-'))):
                stats['improper_strand'] += 1
                continue
            if (strandness == 'reverse') and (not ((read1.iv.strand == '-') and (read2.iv.strand == '+'))):
                stats['improper_strand'] += 1
                continue
            # ignore pairs on different chromosomes
            if read1.iv.chrom != read2.iv.chrom:
                stats['diff_chrom'] += 1
                continue
            pos = junction_positions[read1.iv.chrom]
            if read1.iv.start < pos <= read2.iv.end:
                counts[read1.iv.chrom] += 1
        for key, val in stats.items():
            logger.info('{}: {}'.format(key, val))
    else:
        logger.info('count single-end reads')
        for read in sam:
            # ignore unmapped read
            if not read.aligned:
                continue
            # ignore reads with mapping quality below threshold
            if read.aQual < min_mapping_quality:
                continue
            if (strandness == 'forward') and (read.iv.strand == '-'):
                continue
            if (strandness == 'reverse') and (not ((read.iv.strand == '+'))):
                continue
            pos = junction_positions[read.iv.chrom]
            if read.iv.start < pos <= read.iv.end:
                counts[read.iv.chrom] += 1
    # output counts
    logger.info('count fragments: {}'.format(counts.sum()))
    logger.info('write counts to file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as fout:
        counts.to_csv(fout, sep='\t', header=None, index=True, na_rep='NA')

@command_handler
def count_mature_mirna(args):
    from collections import OrderedDict, defaultdict
    from ioutils import open_file_or_stdin, open_file_or_stdout
    import pysam
    from utils import read_gff

    logger.info('read input GFF file: ' + args.annotation)
    fin = open(args.annotation, 'r')
    # key: precursor_id, value: precursor record
    precursors = OrderedDict()
    # key: precursor_id, value: list of mature records
    matures = defaultdict(list)
    mature_names = []
    # read features from GFF file
    for record in read_gff(fin):
        if record.feature == 'miRNA_primary_transcript':
            precursors[record.attr['ID']] = record
        elif record.feature == 'miRNA':
            matures[record.attr['Derives_from']].append(record)
            mature_names.append(record.attr['Name'])
    fin.close()
    # get locations of mature miRNAs
    # key: precursor_name, key: dict of (mature_name, (start, end))
    mature_locations = defaultdict(dict)
    for precursor_id, precursor in precursors.items():
        precursor_name = precursor.attr['Name']
        for mature in matures[precursor_id]:
            if mature.strand == '+':
                mature_locations[precursor_name][mature.attr['Name']] = (
                    mature.start - precursor.start,
                    mature.end - precursor.start + 1)
            else:
                mature_locations[precursor_name][mature.attr['Name']]  = (
                    precursor.end - mature.end,
                    precursor.end - mature.start + 1)

    logger.info('read input BAM/SAM file: ' + args.input_file)
    sam = pysam.AlignmentFile(args.input_file, "rb")
    counts = defaultdict(int)
    min_mapping_quality = args.min_mapping_quality
    for read in sam:
        if read.is_unmapped:
            continue
        if read.mapping_quality < min_mapping_quality:
            continue
        if read.is_reverse:
            continue
        # find mature miRNA with maximum overlap with the read
        max_overlap = 0
        matched_mature_name = None
        for mature_name, mature_location in mature_locations[read.reference_name].items():
            # get overlap
            overlap = (mature_location[1] - mature_location[0]) \
                + (read.reference_end - read.reference_start) \
                - (max(read.reference_end, mature_location[1]) - min(read.reference_start, mature_location[0]))
            if overlap > max_overlap:
                max_overlap = overlap
                matched_mature_name = mature_name
        if max_overlap <= 0:
            continue
        # count the read
        counts[matched_mature_name] += 1
    
    logger.info('open output file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as f:
        for mature_name in mature_names:
            f.write('{}\t{}\n'.format(mature_name, counts[mature_name]))
         
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
    parser.add_argument('--strandness', '-s', type=str, default='no',
        choices=('forward', 'reverse', 'no'),
        help='forward/reverse: only count reads in reverse strand. no: count reads in both strands')
    parser.add_argument('--output-file', '-o', type=str, default='-', 
        help='output tab-deliminated file. Two columns: gene_id, count')
    
    parser = subparsers.add_parser('count_mature_mirna',
        help='count reads mapped to mature miRNA')
    parser.add_argument('--input-file', '-i', type=str, required=True, 
        help='input BAM/SAM file mapped to miRBase hairpin sequences')
    parser.add_argument('--annotation', '-a', type=str, required=True,
        help='GFF3 file containing mature miRNA locations in precursor miRNA')
    parser.add_argument('--min-mapping-quality', '-q', type=int, default=0,
        help='only count reads with mapping quality greater than this number')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output file')
    
    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('count_reads.' + args.command)

    command_handlers.get(args.command)(args)