#! /usr/bin/env python
import os
import sys
import argparse
import logging
import pickle
import json
from flask import Flask, request, send_from_directory, send_file
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('web_server')

app = Flask(__name__, static_url_path='')

def read_gtf(filename):
    with open(filename, 'r') as fin:
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
            #gene_id = attrs.get('gene_id')
            #if gene_id is None:
            #    raise ValueError('gene_id not found in GTF file at line {}'.format(lineno))
            yield (c, attrs, line, lineno)

def read_gff3(filename):
    with open(filename, 'r') as fin:
        lineno = 0
        for line in fin:
            lineno += 1
            c = line.strip().split('\t')
            if c[0].startswith('#'):
                continue
            attrs = {}
            for a in c[8].split(';'):
                a = a.strip()
                i = a.find('=')
                key = a[:i]
                val = a[(i + 1):].strip('"')
                attrs[key] = val
            yield (c, attrs, line, lineno)

def build_feature_db_from_gtf(filename, source=None):
    genes = {}
    transcripts = {}
    gene_names = {}
    transcript_names = {}
    for c, attrs, line, lineno in read_gtf(filename):
        if c[2] == 'exon':
            gene_id = attrs.get('gene_id')
            if gene_id is None:
                raise ValueError('gene_id attribute not found in GTF file {}:{}'.format(filename, lineno))
            gene = genes.get(gene_id)
            if gene is None:
                # chrom, start, end, strand, source
                gene = [c[0], int(c[3]) - 1, int(c[4]), c[6], c[1]]
                genes[gene_id] = gene
            else:
                gene[1] = min(gene[1], int(c[3]) - 1)
                gene[2] = max(gene[2], int(c[4]))
            
            transcript_id = attrs.get('transcript_id')
            if transcript_id is None:
                raise ValueError('transcript_id attribute not found in GTF file {}:{}'.format(filename, lineno))
            transcript = transcripts.get(transcript_id)
            if transcript is None:
                # chrom, start, end, strand, source
                transcript = [c[0], int(c[3]) - 1, int(c[4]), c[6], c[1]]
                transcripts[transcript_id] = transcript
            else:
                transcript[1] = min(transcript[1], int(c[3]) - 1)
                transcript[2] = max(transcript[2], int(c[4]))
            
            gene_name = attrs.get('gene_name')
            if gene_name is not None:
                gene_names[gene_name] = gene
            transcript_name = attrs.get('transcript_name')
            if transcript_name is not None:
                transcript_names[transcript_name] = transcript
    feature_db = {}
    feature_db.update(genes)
    feature_db.update(transcripts)
    feature_db.update(gene_names)
    feature_db.update(transcript_names)
    return feature_db

def build_feature_db_from_gff3(filename, source=None):
    transcripts = {}
    genes = {}
    aliases = {}
    for c, attrs, line, lineno in read_gff3(filename):
        name = attrs.get('Name')
        gene = attrs.get('gene')
        # try to use gene or Name attribute as gene_id
        if gene is not None:
            gene_id = gene
        elif name is not None:
            gene_id = name
        else:
            raise ValueError('neither Name or gene attribute is not found in GFF3 file {}:{}'.format(filename, lineno))
        # update gene range
        feature = genes.get(gene_id)
        if feature is None:
            feature = [c[0], int(c[3]) - 1, int(c[4]), c[6], c[1]]
            genes[gene_id] = feature
        else:
            feature[1] = min(feature[1], int(c[3]) - 1)
            feature[2] = max(feature[2], int(c[4]))
        # alias
        alias = aliases.get('Alias')
        if alias is not None:
            aliases[alias] = feature
        
        transcript_id = attrs.get('transcript_id')
        feature = transcripts.get(transcript_id)
        if feature is not None:
            # update transcript range
            if feature is None:
                feature = [c[0], int(c[3]) - 1, int(c[4]), c[6], c[1]]
                transcripts[transcript_id] = feature
            else:
                feature[1] = min(feature[1], int(c[3]) - 1)
                feature[2] = max(feature[2], int(c[4]))
  
    feature_db = {}
    #feature_db.update(features)
    feature_db.update(genes)
    feature_db.update(aliases)
    feature_db.update(transcripts)
    return feature_db

def build_feature_db_from_bed(filename, source='BED'):
    feature_db = {}
    with open(filename, 'r') as f:
        for line in f:
            c = line.strip().split('\t')
            feature_db[c[3]] = [c[0], int(c[1]), int(c[2]), c[5], source]
    return feature_db

def build_feature_db_from_genepred(filename, source='GENEPRED'):
    feature_db = {}
    with open(filename, 'r') as f:
        for line in f:
            c = line.strip().split('\t')
            feature_db[c[0]] = [c[1], int(c[3]), int(c[4]), c[2], source]
    return feature_db

def build_feature_db_from_genepredext(filename, source='REFGENE'):
    feature_db = {}
    with open(filename, 'r') as f:
        for line in f:
            c = line.strip().split('\t')
            feature = [c[1], int(c[3]), int(c[4]), c[2], source]
            feature_db[c[1]] = feature
            feature_db[c[11]] = feature
    return feature_db

def build_feature_db_from_refgene(filename, source='REFGENE'):
    feature_db = {}
    with open(filename, 'r') as f:
        for line in f:
            c = line.strip().split('\t')
            feature = [c[2], int(c[4]), int(c[5]), c[3], source]
            feature_db[c[1]] = feature
            feature_db[c[12]] = feature
    return feature_db

def build_feature_db(filename, source=None, format=None):
    if format is None:
        format = filename.split('.')[-1].lower()
    if format == 'bed':
        return build_feature_db_from_bed(filename, source=source)
    elif format == 'gtf':
        return build_feature_db_from_gtf(filename, source=source)
    elif format == 'gff3':
        return build_feature_db_from_gff3(filename, source=source)
    elif format == 'genepred':
        return build_feature_db_from_genepred(filename, source=source)
    elif format == 'genepredext':
        return build_feature_db_from_genepredext(filename, source=source)
    elif format == 'refgene':
        return build_feature_db_from_refgene(filename, source=source)
    elif format == 'pkl':
        return pickle.load(open(filename, 'rb'))
    else:
        raise ValueError('unknown database file format: {}'.format(format))
    
def read_feature_db(filename):
    feature_db = {}
    with open(filename, 'r') as f:
        for line in f:
            c = line.strip().split('\t')
            feature_db[c[0]] = c[1:]
    return feature_db

def merge_feature_db(filenames):
    feature_db = {}
    for filename in filenames:
        with open(filename, 'rb') as f:
            feature_db.update(pickle.load(f))
    return feature_db

@app.route('/locus', methods=['GET'])
def query_locus():
    global feature_db
    genome = request.args.get('genome')
    name = request.args.get('name')
    if (genome is None) or (name is None):
        return '[]'
    db = feature_db.get(genome)
    if db is None:
        return '[]'
    feature = db.get(name)
    if feature is not None:
        #return '{0}\t{1}:{2}-{3}\t{4}'.format(name, feature[0], feature[1], feature[2], feature[4])
        return json.dumps([{
            'gene': name,
            'chromosome': feature[0],
            'start': feature[1],
            'end': feature[2],
            'type': 'gene'
        }])
    else:
        return '[]'
    
@app.route('/igv/<path:filename>')
def serve_static_file(filename):
    return send_from_directory(os.path.join(os.getcwd(), 'igv', 'html'), filename)

@app.route('/bigwig/<dataset>/<filename>')
def serve_bigwig_file(dataset, filename):
    return send_file(os.path.join(os.getcwd(), 'output', dataset, 'bigwig', filename), conditional=True)

@app.route('/genome/hg38/<path:filename>')
def serve_genome_dir(filename):
    return send_from_directory(os.path.join(os.getcwd(), 'genome', 'hg38'), filename, conditional=True)

if __name__ == "__main__":
    #app.run(host='0.0.0.0', debug=True)
    parser = argparse.ArgumentParser('IGV web server')
    parser.add_argument('--database', '-i', type=str, action='append', help='database files')
    parser.add_argument('--host', type=str, default='0.0.0.0')
    parser.add_argument('--port', type=int, default=5000)
    parser.add_argument('--output-file', '-o', type=str, help='feature database file')
    parser.add_argument('--build-database', action='store_true', default=False)
    parser.add_argument('--genome', type=str)
    parser.add_argument('--debug', action='store_true', default=False)
    args = parser.parse_args()

    if args.build_database:
        if (args.genome is None) or (args.database is None):
            print('Error: argument --genome and --database is required for --build-database')
            parser.print_help()
            sys.exit(1)
    feature_db = {}
    if args.database is not None:
        for database in args.database:
            if args.build_database:
                logger.info('build feature database from {}'.format(database))
                feature_db.update(build_feature_db(database))
            else:
                with open(database, 'rb') as f:
                    genome = os.path.splitext(os.path.basename(database))[0]
                    logger.info('load feature database {} from {}'.format(genome, database))
                    feature_db[genome] = pickle.load(f)
                #feature_db.update(read_feature_db(database))
            logger.info('finished building feature database')
    if args.output_file is not None:
        logger.info('dump database to output file: ' + args.output_file)
        with open(args.output_file, 'wb') as f:
            pickle.dump(feature_db, f)
    if not args.build_database:
        app.run(host=args.host, port=args.port, debug=args.debug, threaded=True)