#! /usr/bin/env python
from __future__ import print_function
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

def read_gtf(filename):
    from ioutils import open_file_or_stdin

    with open_file_or_stdin(filename) as fin:
        lineno = 0
        for line in fin:
            lineno += 1
            c = line.strip().split('\t')
            if c[0].startswith('#'):
                continue
            attrs = {}
            for a in c[8].split(';')[:-1]:
                a = a.strip()
                i = a.find(' ')
                key = a[:i]
                val = a[(i + 1):].strip('"')
                attrs[key] = val
            gene_id = attrs.get('gene_id')
            if gene_id is None:
                raise ValueError('gene_id not found in GTF file at line {}'.format(lineno))
            yield (c, attrs, line)

@command_handler
def extract_gene(args):
    from ioutils import open_file_or_stdout

    feature = args.feature
    genes = {}
    logger.info('read GTF file: ' + args.input_file)
    for c, attrs, line in read_gtf(args.input_file):
        if c[2] != feature:
            continue
        gene_id = attrs.get('gene_id')
        gene = genes.get(gene_id)
        if gene is None:
            gene = [c[0], int(c[3]) - 1, int(c[4]), gene_id, 0, c[6]]
            genes[gene_id] = gene
        else:
            gene[1] = min(gene[1], int(c[3]) - 1)
            gene[2] = max(gene[2], int(c[4]))
    
    logger.info('number of genes: {}'.format(len(genes)))
    logger.info('write BED file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as fout:
        for gene_id, gene in genes.items():
            fout.write('\t'.join(map(str, gene)))
            fout.write('\n')


@command_handler
def fix_gtf(args):
    from ioutils import open_file_or_stdout
    from collections import defaultdict
    # strand of exons grouped by transcript_id
    strands = defaultdict(list)

    feature = args.feature
    lines = []
    logger.info('read GTF file: ' + args.input_file)
    for c, attrs, line in read_gtf(args.input_file):
        if c[2] in ('transcript', 'exon'):
            transcript_id = attrs.get('transcript_id')
            if transcript_id is None:
                raise ValueError('transcript_id not found in GTF file at line {}'.format(lineno))
            lines.append((transcript_id, line))
        else:
            transcript_id = attrs.get('transcript_id')
            if transcript_id is None:
                raise ValueError('transcript_id not found in GTF file at line {}'.format(lineno))
            strands[transcript_id].append(c[6])
    invalid_transcripts = set()
    for transcript_id, strands_tx in strands.items():
        strands_tx = set(strands_tx)
        # remove transcripts without strand information
        if '.' in strands_tx:
            invalid_transcripts.add(transcript_id)
        # remove transcripts with exon on different strands
        elif len(strands_tx) != 1:
            invalid_transcripts.add(transcript_id)
    
    logger.info('number of transcripts: {}'.format(len(strands)))
    logger.info('number of invalid transcripts: {}'.format(len(invalid_transcripts)))
    logger.info('write GTF file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as fout:
        for transcript_id, line in lines:
            if transcript_id not in invalid_transcripts:
                fout.write(line)


@command_handler
def transcript_counts(args):
    import pysam
    import numpy as np
    from ioutils import open_file_or_stdout
    from collections import defaultdict

    logger.info('read input transcript BAM file: ' + args.input_file)
    sam = pysam.AlignmentFile(args.input_file, "rb")
    counts = defaultdict(int)
    for read in sam:
        counts[read.reference_name] += 1
    
    logger.info('create output file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as fout:
        for key, val in counts.items():
            fout.write('{}\t{}\n'.format(key, val))

@command_handler
def extract_longest_transcript(args):
    from ioutils import open_file_or_stdin, open_file_or_stdout
    from collections import defaultdict
    from functools import partial

    feature = args.feature
    genes = defaultdict(partial(defaultdict, int))
    lines = []
    logger.info('read gtf file: ' + args.input_file)
    with open_file_or_stdin(args.input_file) as fin:
        lineno = 0
        for line in fin:
            lineno += 1
            c = line.strip().split('\t')
            if c[0].startswith('#'):
                continue
            if c[2] != feature:
                lines.append(('#other#', line))
                continue
            attrs = {}
            for a in c[8].split(';')[:-1]:
                a = a.strip()
                i = a.find(' ')
                key = a[:i]
                val = a[(i + 1):].strip('"')
                attrs[key] = val
            transcript_id = attrs.get('transcript_id')
            if transcript_id is None:
                raise ValueError('transcript_id not found in GTF file at line {}'.format(lineno))
            gene_id = attrs.get('gene_id')
            if gene_id is None:
                raise ValueError('gene_id not found in GTF file at line {}'.format(lineno))
            lines.append((transcript_id, line))
            genes[gene_id][transcript_id] += int(c[4]) - int(c[3]) + 1
    kept_transcripts = set()
    kept_transcripts.add('#other#')
    for gene_id, gene in genes.items():
        max_length = 0
        max_transcript = None
        for transcript_id, length in gene.items():
            if length > max_length:
                max_length = length
                max_transcript = transcript_id
        kept_transcripts.add(transcript_id)

    logger.info('number of genes: {}'.format(len(genes)))
    logger.info('number of transcripts: {}'.format(sum(map(len, genes.values()))))
    logger.info('number of longest transcripts: {}'.format(len(kept_transcripts) - 1))
    logger.info('write output gtf file: ' + args.output_file)
    with open_file_or_stdout(args.output_file) as fout:
        for transcript_id, line in lines:
            if transcript_id in kept_transcripts:
                fout.write(line)

@command_handler
def gtf_to_transcript_table(args):
    from ioutils import open_file_or_stdin, open_file_or_stdout
    from collections import OrderedDict

    feature = args.feature
    default_transcript_type = args.transcript_type
    default_gene_type = args.gene_type

    fout = open_file_or_stdout(args.output_file)
    with open_file_or_stdin(args.input_file) as fin:
        transcripts = OrderedDict()
        for line in fin:
            c = line.strip().split('\t')
            if c[0].startswith('#'):
                continue
            if c[2] != feature:
                continue
            attrs = {}
            for a in c[8].split(';')[:-1]:
                a = a.strip()
                i = a.find(' ')
                key = a[:i]
                val = a[(i + 1):].strip('"')
                attrs[key] = val
            if 'transcript_name' not in attrs:
                attrs['transcript_name'] = attrs['transcript_id']
            if 'gene_name' not in attrs:
                attrs['gene_name'] = attrs['gene_id']
            if default_transcript_type is not None:
                attrs['transcript_type'] = default_transcript_type
            else:
                if 'transcript_type' not in attrs:
                    attrs['transcript_type'] = 'unknown'
            if default_gene_type is not None:
                attrs['gene_type'] = default_gene_type
            else:
                if 'gene_type' not in attrs:
                    attrs['gene_type'] = 'unknown'
            exon = [c[0], int(c[3]) - 1, int(c[4]), attrs['gene_id'], 0, c[6],
                attrs['gene_id'], attrs['transcript_id'], 
                attrs['gene_name'], attrs['transcript_name'],
                attrs['gene_type'], attrs['transcript_type'], c[1]]
            transcript = transcripts.get(attrs['transcript_id'])
            if transcript is None:
                transcripts[attrs['transcript_id']] = exon
            else:
                if c[2] == 'exon':
                    transcript[1] = min(transcript[1], exon[1])
                    transcript[2] = max(transcript[2], exon[2])
        header = ['chrom', 'start', 'end', 'name', 'score', 'strand',
            'gene_id', 'transcript_id', 
            'gene_name', 'transcript_name',
            'gene_type', 'transcript_type', 'source'
        ]
        print('\t'.join(header), file=fout)
        for transcript in transcripts.values():
            print('\t'.join(str(a) for a in transcript), file=fout)
    fout.close()

@command_handler
def extract_circrna_junction(args):
    from Bio import SeqIO
    from ioutils import open_file_or_stdin, open_file_or_stdout

    anchor_size = args.anchor_size
    logger.info('read sequence file: ' + args.input_file)
    logger.info('create output file: ' + args.output_file)
    fout = open_file_or_stdout(args.output_file)
    with open_file_or_stdin(args.input_file) as fin:
        for record in SeqIO.parse(fin, 'fasta'):
            seq = str(record.seq)
            if len(seq) < args.min_length:
                continue
            s = min(len(seq), anchor_size)
            seq_id = record.id.split('|')[0]
            fout.write('>{}\n'.format(seq_id))
            fout.write(seq[-s:] + seq[:s])
            fout.write('\n')
    fout.close()

@command_handler
def filter_circrna_reads(args):
    import pysam
    import numpy as np
    from ioutils import open_file_or_stdout, open_file_or_stdin
    from collections import defaultdict
    from copy import deepcopy

    logger.info('read input SAM file: ' + args.input_file)
    fin = open_file_or_stdin(args.input_file)
    sam_in = pysam.AlignmentFile(fin, "r")
    if sam_in.header is None:
        raise ValueError('requires SAM header to get junction positions')
    # get junction positions (middle of the sequences)
    junction_positions = {}
    for sq in sam_in.header['SQ']:
        junction_positions[sq['SN']] = sq['LN']//2

    logger.info('create output SAM file: ' + args.output_file)
    fout = open_file_or_stdout(args.output_file)
    sam_out = pysam.AlignmentFile(fout, 'w', template=sam_in)

    sam_filtered = None
    if args.filtered_file is not None:
        logger.info('create filtered SAM file: ' + args.filtered_file)
        sam_filtered = pysam.AlignmentFile(args.filtered_file, 'w', template=sam_in)

    for read in sam_in:
        filtered = False
        if read.is_unmapped:
            filtered = True
        elif read.is_reverse:
            filtered = True
        else:
            pos = junction_positions[read.reference_name]
            if not (read.reference_start < pos <= read.reference_end):
                filterd = True
        if not filtered:
            sam_out.write(read)
        elif sam_filtered is not None:
            sam_filtered.write(read)
    
    fin.close()
    fout.close()
    if sam_filtered is not None:
        sam_filtered.close()

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Preprocessing module')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('transcript_counts', 
        help='count reads that belongs to a transcript')
    parser.add_argument('--input-file', '-i', type=str, required=True,
        help='input transcript BAM file')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output transcript counts file')
    
    parser = subparsers.add_parser('gtf_to_transcript_table')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='input GTF file')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output table file')
    parser.add_argument('--feature', type=str,
        help='feature to use in input GTF file (Column 3)')
    parser.add_argument('--gene-type', type=str,
        help='gene type to set if "gene_type" attribute is not available in GTF file')
    parser.add_argument('--transcript-type', type=str, default='unknown',
        help='gene type to set if "transcript_type" attribute is not available in GTF file')
    
    parser = subparsers.add_parser('extract_longest_transcript')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='input GTF file')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output table file')
    parser.add_argument('--feature', type=str, default='exon',
        help='feature to use in input GTF file (Column 3)')
    
    parser = subparsers.add_parser('fix_gtf',
        help='fix problems in a GTF file by removing invalid transcripts')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='input GTF file')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output table file')
    parser.add_argument('--feature', type=str, default='exon',
        help='feature to use in input GTF file (Column 3)')
    
    parser = subparsers.add_parser('extract_gene',
        help='extract gene regions from a GTF file')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='input GTF file')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output BED file')
    parser.add_argument('--feature', type=str, default='exon',
        help='feature to use in input GTF file (Column 3)')
    
    parser = subparsers.add_parser('extract_circrna_junction',
        help='extract circular RNA junction sequences from spliced sequences')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='spliced sequence file')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='junction sequence file')
    parser.add_argument('--anchor-size', '-s', type=int, default=50,
        help='anchor length from both ends')
    parser.add_argument('--min-length', type=int, default=50,
        help='minimal length to keep a spliced sequence')
    
    parser = subparsers.add_parser('filter_circrna_reads',
        help='filter out reads not crossing circRNA junctions')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='input SAM file')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output SAM file')
    parser.add_argument('--filtered-file', '-u', type=str,
        help='write filtered SAM records to file')

    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('preprocess.' + args.command)
    command_handlers.get(args.command)(args)