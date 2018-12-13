#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

@command_handler
def sequential_mapping(args):
    from jinja2 import Template, Environment
    from ioutils import open_file_or_stdout

    rna_types = []
    if len(args.rna_types) > 0:
        rna_types = args.rna_types.split(',')

    logger.info('load template: ' + args.template)
    env = Environment(lstrip_blocks=True, trim_blocks=True)
    with open(args.template, 'r') as f:
        template = env.from_string(f.read())
    with open_file_or_stdout(args.output_file) as f:
        f.write(template.render(rna_types=rna_types, aligner=args.aligner))

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Count reads in BAM files')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('sequential_mapping', 
        help='generate Snakefile for sequential mapping')
    parser.add_argument('--rna-types', type=str, required=True,
        help='comma-separated list of rna types')
    parser.add_argument('--aligner', type=str, default='bowtie2',
        help='aligner to use')
    parser.add_argument('--template', '-t', type=str, default='templates/sequential_mapping.snakemake',
        help='template for snakefile')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output file')
    
    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('generate_snakemake.' + args.command)

    command_handlers.get(args.command)(args)