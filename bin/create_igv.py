#! /usr/bin/env python
import argparse, sys, os, errno
import logging
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

    


if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Create IGV')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('create_reference', 
        help='Create IGV reference genome from FASTA file')
    parser.add_argument('--genome', type=str, help='genome name')
    parser.add_argument('--name', type=str, help='reference name')
    parser.add_argument('--fasta', type=str, required=True)
    parser.add_argument('--fasta-index', type=str, required=True)
    parser.add_argument('--gene-name', type=str,
        help='Tab-separated file with 2 columns: gene_id, gene_name')
    parser.add_argument('--gene-id', type=str,
        help='Only use gene ids in provided list file')
    parser.add_argument('--output-dir', '-o', type=str,
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
