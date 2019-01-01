#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
    
def create_igv(args):
    from jinja2 import Template, Environment
    from ioutils import open_file_or_stdout

    env = Environment(lstrip_blocks=True, trim_blocks=True)
    with open(args.input_file, 'r') as f:
        template = env.from_string(f.read())
    data = {
        'genome': args.genome,
        'locus': args.locus,
        'reference': None
    }

    if args.genome_dir is not None:
        data['reference'] = {
            'id': args.genome,
            'fastaURL': os.path.join(args.base_url, args.genome_dir, 'fasta', 'genome.fa'),
            'indexURL': os.path.join(args.base_url, args.genome_dir, 'fasta', 'genome.fa.fai'),
            'cytobandURL': os.path.join(args.base_url, args.genome_dir, 'igv', 'cytoBandIdeo.txt'),
            'geneURL': os.path.join(args.base_url, args.genome_dir, 'igv', 'refGene.txt')
        }
    
    data['tracks'] = [
        {'type': 'annotation', 
        'name': 'refGene', 
        'format': 'refgene', 
        'url': os.path.join(args.base_url, args.genome_dir, 'igv', 'refGene.txt')
        }
    ]
    file_types = {
        'bigWig': {'type': 'wig', 'format': 'bigwig'},
        'gff3': {'type': 'annotation', 'format': 'gff3'},
        'bed': {'type': 'annotation', 'format': 'bed'},
        'genePredExt': {'type': 'annotation', 'format': 'genepredext'},
        'genePred': {'type': 'annotation', 'format': 'genepred'}
    }
    #data['tracks'] = []
    if args.track:
        for track in args.track:
            track_meta = {}
            for file_type in file_types.keys():
                if track.endswith(file_type):
                    track_meta.update(file_types[file_type])
            track_meta['name'] = os.path.basename(track)
            track_meta['url'] = os.path.join(args.base_url, track)
            if os.path.exists(track + '.idx'):
                track_meta['indexURL'] = os.path.join(args.base_url, track + '.idx')
            data['tracks'].append(track_meta)

    with open_file_or_stdout(args.output_file) as f:
        f.write(template.render(**data))
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create an IGV web page')
    parser.add_argument('--input-file', '-i', type=str, required=True,
        help='template file for IGV files')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output file')
    parser.add_argument('--genome-dir', '-r', type=str, required=True,
        help='reference genome directory from http://s3.amazonaws.com/igv.broadinstitute.org/genomes')
    parser.add_argument('--track', type=str, action='append',
        help='track file in any supported format in https://software.broadinstitute.org/software/igv/RecommendedFileFormats')
    parser.add_argument('--genome', '-g', type=str, required=True)
    parser.add_argument('--locus', '-l', type=str, default='',
        help='genomic locus, e.g. chr1:1000000-1000100')
    parser.add_argument('--base-url', type=str, default='.')

    args = parser.parse_args()
    create_igv(args)
