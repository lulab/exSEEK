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
def gff3_to_transcript_table(args):
    from ioutils import open_file_or_stdin, open_file_or_stdout
    from collections import OrderedDict

    fout = open_file_or_stdout(args.output_file)
    with open_file_or_stdin(args.input_file) as fin:
        transcripts = OrderedDict()
        for line in fin:
            if line.startswith('#'):
                continue
            c = line.strip().split('\t')
            attrs = {}
            for a in c[8].split(';'):
                key, value = a.split('=')
                attrs[key] = value
            


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
            # mitranscriptome
            if c[1] == 'mitranscriptome':
                if attrs['tcat'] == 'lncrna':
                    attrs['gene_type'] = 'lncRNA'
                    attrs['transcript_type'] = 'lncRNA'
                elif attrs['tcat'] == 'tucp':
                    attrs['gene_type'] = 'tucpRNA'
                    attrs['transcript_type'] = 'tucpRNA'
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

@command_handler
def chrom_sizes(args):
    from Bio import SeqIO
    from ioutils import open_file_or_stdin, open_file_or_stdout

    fout = open_file_or_stdout(args.output_file)
    with open_file_or_stdin(args.input_file) as fin:
        for record in SeqIO.parse(fin, 'fasta'):
            fout.write('{}\t{}\n'.format(record.id, len(record.seq)))

@command_handler
def extract_feature_sequence(args):
    from pyfaidx import Fasta
    from Bio.Seq import Seq
    from ioutils import open_file_or_stdout

    fout = open_file_or_stdout(args.output_file)
    fastas = {}
    with open(args.input_file, 'r') as fin:
        for lineno, line in enumerate(fin):
            feature = line.split('\t')[0]
            gene_id, gene_type, gene_name, domain_id, transcript_id, start, end = feature.split('|')
            start = int(start)
            end = int(end)
            if gene_type == 'genomic':
                gene_type = 'genome'
            if gene_type not in fastas:
                fastas[gene_type] = Fasta(os.path.join(args.genome_dir, 'fasta', gene_type + '.fa'))
            if gene_type == 'genome':
                chrom, gstart, gend, strand = gene_id.split('_')
                gstart = int(gstart)
                gend = int(gend)
                seq = fastas[gene_type][chrom][gstart:gend].seq
                if strand == '-':
                    seq = str(Seq(seq).reverse_complement())
            else:
                seq = fastas[gene_type][transcript_id][start:end].seq
            seq = seq.upper()
            fout.write('>{}\n'.format(feature))
            fout.write(seq)
            fout.write('\n')
    fout.close()

@command_handler
def calculate_star_parameters(args):
    from math import log2

    genome_length = 0
    n_seqs = 0
    with open(args.input_file, 'r') as f:
        for line in f:
            genome_length += int(line.split('\t')[1])
            n_seqs += 1
    if args.parameter == 'genomeSAindexNbases':
        print(min(14, int(log2(genome_length)//2) - 1))
    elif args.parameter == 'genomeChrBinNbits':
        print(min(18, int(log2(genome_length/n_seqs))))

@command_handler
def get_mature_mirna_location(args):
    import gzip
    from Bio import SeqIO
    from collections import defaultdict
    from ioutils import open_file_or_stdout

    logger.info('read hairpin sequences: ' + args.hairpin)
    hairpin_seqs = {}
    with open(args.hairpin, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if (args.species is None) or (args.species == record.id.split('-')[0]):
                hairpin_seqs[record.id] = str(record.seq)
    
    logger.info('read mature sequences: ' + args.mature)
    fout = open_file_or_stdout(args.output_file)
    mature_positions = defaultdict(dict)
    with open(args.mature, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            if not ((args.species is None) or (args.species == record.id.split('-')[0])):
                continue
            if record.id.endswith('-5p') or record.id.endswith('-3p'):
                end_type = record.id.split('-')[-1]
                hairpin_id = record.id[:-3].lower()
                hairpin_seq = hairpin_seqs.get(hairpin_id)
                if hairpin_seq is None:
                    logger.warn('hairpin sequence id {} is not found in hairpin file'.format(hairpin_id))
                    continue
                mature_id = record.id
                mature_seq = str(record.seq)
                pos = hairpin_seq.find(mature_seq)
                if pos < 0:
                    logger.warn('mature sequence {} is not found in hairpin file'.format(mature_id))
                    continue
                else:
                    mature_positions[hairpin_id][end_type] = (pos, pos + len(mature_seq))

    for hairpin_id, hairpin_seq in hairpin_seqs.items():
        mature_position = mature_positions.get(hairpin_id)
        if mature_position is None:
            mature_position = {}
        fout.write('>{}\n'.format(hairpin_id))
        offset = 0
        if '5p' in mature_position:
            start, end = mature_position['5p']
            fout.write('{0}\x1B[31;1m{1}\x1B[0m'.format(hairpin_seq[offset:start], hairpin_seq[start:end]))
            offset = end
        if '3p' in mature_position:
            start, end = mature_position['3p']
            fout.write('{0}\x1B[32;1m{1}\x1B[0m'.format(hairpin_seq[offset:start], hairpin_seq[start:end]))
            offset = end
        fout.write('{}\n'.format(hairpin_seq[offset:]))
    fout.close()

@command_handler
def calculate_gene_length(args):
    import HTSeq
    from collections import defaultdict
    from functools import partial
    import numpy as np
    from ioutils import open_file_or_stdin
    from tqdm import tqdm

    fin = open_file_or_stdin(args.input_file)
    gff = HTSeq.GFF_Reader(fin)
    exons = defaultdict(partial(defaultdict, int))
    for feature in tqdm(gff, unit='feature'):
        if feature.type == 'exon':
            exons[feature.attr['gene_id']][feature.attr['transcript_id']] += feature.iv.length
        
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
    
    parser = subparsers.add_parser('chrom_sizes',
        help='create chrom sizes file from FASTA file')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='input FASTA')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output chrom sizes file')
    
    parser = subparsers.add_parser('extract_feature_sequence',
        help='extract sequences using feature names')
    parser.add_argument('--input-file', '-i', type=str, required=True,
        help='feature file with feature names in the first column')
    parser.add_argument('--genome-dir', '-g', type=str, required=True,
        help='genome directory where fasta/${rna_type}.fa contains sequences')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output file in FASTA format')
    
    parser = subparsers.add_parser('calculate_star_parameters',
        help='calculate --genomeSAindexNbases and --genomeChrBinNbits for STAR')
    parser.add_argument('--input-file', '-i', type=str, required=True,
        help='FASTA index (.fai) file')
    parser.add_argument('--parameter', '-p', type=str, required=True,
        choices=('genomeSAindexNbases', 'genomeChrBinNbits'),
        help='parameter to calculate')

    parser = subparsers.add_parser('get_mature_mirna_location',
        help='get mature miRNA location in precursor miRNA in miRBase')
    parser.add_argument('--mature', type=str, required=True,
        help='FASTA file of mature sequences')
    parser.add_argument('--hairpin', type=str, required=True,
        help='FASTA file of hairpin sequences')
    parser.add_argument('--output-file', '-o', type=str, default='-')
    parser.add_argument('--species', type=str)

    parser = subparsers.add_parser('calculate_gene_length',
        help='calculate effective gene length')
    parser.add_argument('--input-file', '-i', type=str, default='-',
        help='GTF/GFF file')
    parser.add_argument('--method', '-m', type=str,
        choices=('isoform_length_mean', 'isoform_length_median', 'isoform_length_max', 'merged_exon_length'))
    parser.add_argument('--output-file', '-o', type=str,
        help='output tab-deliminated file with two columns: gene_id, length')
    
    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('preprocess.' + args.command)
    command_handlers.get(args.command)(args)