#! /usr/bin/env python
import argparse, sys, os, errno
import logging
import shutil
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

def color_palette(n_colors):
    # list(map(matplotlib.colors.to_hex, sns.color_palette('deep', 10)))
    colors = ['#4c72b0', '#dd8452', '#55a868', '#c44e52', '#8172b3',
     '#937860', '#da8bc3', '#8c8c8c', '#ccb974', '#64b5cd']
    palette = [colors[i%len(colors)] for i in range(n_colors)]
    return palette

@command_handler
def generate_config(args):
    import pandas as pd
    import yaml
    from jinja2 import Template, Environment
    from collections import defaultdict, OrderedDict

    sample_classes = pd.read_table(args.sample_classes, sep='\t', index_col=0, dtype='str').iloc[:, 0]
    # set track colors
    unique_classes = sample_classes.unique()
    config = {}
    config.update(yaml.load(open(args.reference, 'r')))
    for key in ('reference', 'genome'):
        if key not in config:
            raise KeyError('key "{}" is not defined in reference file: {}'.format(key, args.reference))

    config['locus'] = args.locus
    # define groups
    config['groups'] = {}
    colors = {}
    
    for sample_class, color in zip(unique_classes, color_palette(len(unique_classes))):
        colors[sample_class] = color
        config['groups'][sample_class] = dict(color=color, show=True)
    # add tracks
    order = 1
    if 'tracks' not in config:
        config['tracks'] = {}
    else:
        for name, track in config['tracks'].items():
            track['order'] = order
            order += 1
    
    strands = ['+', '-']
    if args.strand is not None:
        strands = [args.strand]
    for sample_class in unique_classes:
        for sample_id in sample_classes[sample_classes == sample_class].index.values[:args.max_samples_per_class]:
            for strand in strands:
                config['tracks']['{0}({1})'.format(sample_id, strand)] = dict(
                    name='{0}({1})'.format(sample_id, strand),
                    sample_id=sample_id,
                    url=args.bigwig_pattern.format(sample_id=sample_id, strand=strand),
                    type='wig',
                    format='bigwig',
                    height=25,
                    group=sample_class,
                    strand=strand,
                    order=order,
                    color=colors[sample_class],
                    autoscale=True,
                    show=True
                )
                order += 1
    # sort tracks
    config['tracks'] = dict(sorted(config['tracks'].items(), key=lambda x: x[1]['order']))
    # add base_url
    for key in ('fastaURL', 'indexURL', 'cytobandURL'):
        if key in config['reference']:
            config['reference'][key] = args.base_url + '/' + config['reference'][key]
    if 'tracks' in config:
        for name, track in config['tracks'].items():
            for key in ('url', 'indexURL'):
                if key in track:
                    track[key] = args.base_url + '/' + track[key]
    
    with open(args.output_file, 'w') as fout:
        yaml.dump(config, fout, default_flow_style=False)

@command_handler
def render(args):
    from jinja2 import Template, Environment, StrictUndefined
    from ioutils import open_file_or_stdout
    from collections import defaultdict
    import yaml

    env = Environment(lstrip_blocks=True, trim_blocks=True, undefined=StrictUndefined)
    with open(args.input_file, 'r') as f:
        template = env.from_string(f.read())

    config = yaml.load(open(args.config, 'r'))
    config['tracks'] = dict(sorted(config['tracks'].items(), key=lambda x: x[1]['order']))
    with open_file_or_stdout(args.output_file) as f:
        f.write(template.render(**config))

@command_handler
def create_reference(args):
    import yaml
    import pandas as pd
    import numpy as np
    from pyfaidx import Fasta
    from Bio.Seq import Seq
    import subprocess

    logger.info('read fasta')
    fasta = Fasta(args.fasta)

    logger.info('read fasta index')
    fasta_index = pd.read_table(args.fasta + '.fai', sep='\t', header=None, dtype=str).iloc[:, [0, 1]].copy()
    fasta_index.columns = ['name', 'length']
    fasta_index.index = fasta_index['name']
    fasta_index.iloc[:, 1] = fasta_index.iloc[:, 1].astype('int')
    # rename gene ids to gene names
    gene_names = None
    if args.gene_names is not None:
        logger.info('convert gene ids to gene names')
        gene_names = pd.read_table(args.gene_names, sep='\t', header=None, dtype=str, index_col=0).iloc[:, 0]
    gene_ids = None
    if args.gene_ids is not None:
        logger.info('read gene ids')
        gene_ids = pd.read_table(args.gene_ids, dtype=str, header=None).iloc[:, 0]

    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir)
    config = {}
    config['genome'] = args.genome
    config['reference'] = {}
    if args.name is not None:
        config['reference']['id'] = args.genome
        if args.name is not None:
            config['reference']['name'] = args.name
        else:
            config['reference']['name'] = args.genome
        config['reference']['fastaURL'] = os.path.join(args.output_dir, 'reference.fa')
        config['reference']['indexURL'] = os.path.join(args.output_dir, 'reference.fa.fai')
        config['reference']['cytobandURL'] = os.path.join(args.output_dir, 'cytoband.txt')
    config['tracks'] = {
        args.genome: {
            'name': args.genome,
            'type': 'annotation',
            'format': 'genepred',
            'url': os.path.join(args.output_dir, 'annotation.genePred'),
            'displayMode': 'EXPANDED',
            'searchable': True,
            'visibilityWindow':  300000000,
            'height': 100,
            'logScale': false,
            'show': True
        }
    }
    if (gene_ids is None) and (gene_names is None):
        logger.info('copy fasta file')
        shutil.copyfile(args.fasta, config['reference']['fastaURL'])
        logger.info('copy fasta index')
        shutil.copyfile(args.fasta + '.fai', config['reference']['indexURL'])
    else:
        logger.info('write subset of reference fasta')
        if gene_ids is None:
            gene_ids = list(fasta.keys())
        with open(config['reference']['fastaURL'], 'w') as fout:
            for gene_id in gene_ids:
                record = fasta[gene_id]
                if record is None:
                    raise ValueError('gene id {} is not found in fasta file'.format(gene_id))
                seq = record[:].seq
                if gene_names is not None:
                    gene_id = gene_names[gene_id]
                fout.write('>{}\n'.format(gene_id))
                fout.write(seq)
                fout.write('\n')
        logger.info('generate fasta index')
        subprocess.check_call(['samtools', 'faidx', config['reference']['fastaURL']])
        fasta_index = fasta_index.loc[gene_ids]
        if gene_names is not None:
            fasta_index['name'] = gene_names[fasta_index['name'].values]

    cytoband = pd.DataFrame(index=np.arange(fasta_index.shape[0]), columns=np.arange(5))
    cytoband.iloc[:, 0] = fasta_index['name'].values
    cytoband.iloc[:, 1] = np.zeros(fasta_index.shape[0], dtype='int')
    cytoband.iloc[:, 2] = fasta_index['length'].values
    cytoband.iloc[:, 3] = ['p%d.1'%i for i in range(fasta_index.shape[0])]
    cytoband.iloc[:, 4] = np.full(fasta_index.shape[0], 'gneg')
    logger.info('write cytoband file')
    cytoband.to_csv(config['reference']['cytobandURL'], sep='\t', header=False, index=False)

    bed = pd.DataFrame(index=np.arange(fasta_index.shape[0]), columns=np.arange(6))
    bed.iloc[:, 0] = fasta_index['name'].values
    bed.iloc[:, 1] = np.zeros(fasta_index.shape[0], dtype='int')
    bed.iloc[:, 2] = fasta_index['length'].values
    bed.iloc[:, 3] = fasta_index['name'].values
    bed.iloc[:, 4] = np.full(fasta_index.shape[0], '0')
    bed.iloc[:, 5] = np.full(fasta_index.shape[0], '+')
    logger.info('write annotation bed file')
    bed_file = os.path.join(args.output_dir, 'annotation.bed')
    bed.to_csv(bed_file, sep='\t', header=False, index=False)
    logger.info('write annotation genePred file')
    subprocess.check_call(['bedToGenePred', bed_file, config['tracks'][args.genome]['url']])

    with open(os.path.join(args.output_dir, 'config.yaml'), 'w') as fout:
        config['reference']['fastaURL'] = config['reference']['fastaURL']
        yaml.dump(config, fout, default_flow_style=False)

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Create IGV')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('create_reference', 
        help='Create IGV reference genome from FASTA file')
    parser.add_argument('--genome', type=str, help='genome name', required=True)
    parser.add_argument('--name', type=str, help='reference name')
    parser.add_argument('--fasta', type=str, required=True)
    parser.add_argument('--gene-names', type=str,
        help='Tab-separated file with 2 columns: gene_id, gene_name')
    parser.add_argument('--gene-ids', type=str,
        help='Only use gene ids in provided list file')
    parser.add_argument('--output-dir', '-o', type=str, required=True,
        help='output directory')

    parser = subparsers.add_parser('render', 
        help='Render HTML from config file')
    parser.add_argument('--input-file', '-i', type=str, required=True,
        help='template file for IGV files')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output file')
    parser.add_argument('--config', '-c', type=str, help='track list in YAML format')

    parser = subparsers.add_parser('generate_config', 
        help='Generate config from track files')
    parser.add_argument('--sample-classes', type=str, required=True,
        help='sample classes file')
    parser.add_argument('--reference', type=str, required=True,
        help='config file for reference genome')
    parser.add_argument('--max-samples-per-class', type=int, default=10)
    parser.add_argument('--bigwig-pattern', type=str, required=True,
        help='format string with one variables; sample_id, strand')
    parser.add_argument('--locus', '-l', type=str, default='',
        help='genomic locus, e.g. chr1:1000000-1000100')
    parser.add_argument('--base-url',  type=str, default='.')
    parser.add_argument('--strand', '-s', type=str)
    parser.add_argument('--output-file', '-o', type=str, required=True)

    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('create_igv.' + args.command)

    try:
        command_handlers.get(args.command)(args)
    except BrokenPipeError:
        pass
    except KeyboardInterrupt:
        pass
