#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

def get_aligner_index(aligner, rna_type, genome_dir=''):
    if rna_type == 'genome':
        return genome_dir + '/genome_index/{aligner}/genome'.format(aligner=aligner)

    if aligner == 'bowtie2':
        return genome_dir + '/rsem_index/{aligner}/{rna_type}'.format(
            aligner=aligner, rna_type=rna_type)
    elif aligner == 'star':
        return genome_dir + '/rsem_index/{aligner}/{rna_type}'.format(
            aligner=aligner, rna_type=rna_type)
    else:
        raise ValueError('unknown aligner: {}'.format(aligner))

@command_handler
def sequential_mapping(args):
    import string
    from ioutils import open_file_or_stdout

    rna_types = args.rna_types.split(',')

    if args.aligner == 'bowtie2':
        rule_template = string.Template("""
rule map_${rna_type}:
    input:
        reads='${reads}',
        index=genome_dir + '${index}.1.bt2'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/${rna_type}.fa.gz',
        bam='{output_dir}/${bam_type}/{sample_id}/${rna_type}.bam'
    params:
        index=genome_dir + '${index}'
    threads: 
        config['threads_mapping']
    shell:
        '''pigz -d -c {input.reads} \
        | bowtie2 -f -p {threads} --norc --sensitive --no-unal \\
            --un-gz {output.unmapped} -x {params.index} - -S - \\
        | samtools view -b -o {output.bam}
        '''
""")
    elif args.aligner == 'star':
        rule_template = string.Template("""
rule map_${rna_type}:
    input:
        reads='${reads}',
        index=genome_dir + '${index}/SA'
    output:
        unmapped='{output_dir}/unmapped/{sample_id}/${rna_type}.fa',
        bam='{output_dir}/${bam_type}/{sample_id}/${rna_type}.bam',
        star_bam='{output_dir}/star_small/{sample_id}/${rna_type}/Aligned.out.bam'
    params:
        index=genome_dir + '${index}'
        output_prefix='{output_dir}/star_small/{sample_id}/${rna_type}',
        star_unmapped='{output_dir}/star_small/{sample_id}/${rna_type}/Unmapped.out'
    threads:
        config['threads_mapping']
    shell:
        '''STAR --genomeDir {params.index} \
            --readFilesIn {input.reads} \
            --readFilesCommand pigz -d -c \
            --runThreadN {threads} \
            --outReadsUnmapped Fastx \
            --outputFileNamePrefix {params.output_prefix} \
            --outputSAMtype BAM Unsorted
            ln -f {output.star_bam} {output.bam}
            pigz -c {params.star_unmapped} > {output.unmapped}
            rm -f {params.star_unmapped}
        '''
""")
    else:
        raise ValueError('unknown aligner: {}'.format(wildcards.aligner))
    reads = '{output_dir}/unmapped/{sample_id}/clean.fa.gz'
    content = ''
    # sequential mapping
    for i, rna_type in enumerate(rna_types):
        index = get_aligner_index(aligner=args.aligner, rna_type=rna_type)
        if i > 0:
            reads = '{{output_dir}}/unmapped/{{sample_id}}/{rna_type}.fa.gz'.format(rna_type=rna_types[i - 1])
        content = content + rule_template.substitute(reads=reads, rna_type=rna_type, index=index, bam_type='tbam')
    # map to other genomic locations
    if len(rna_types) > 0:
        index = get_aligner_index(aligner=args.aligner, rna_type='genome')
        content = content + rule_template.substitute(reads=reads, rna_type='other', index=index, bam_type='gbam')
    with open_file_or_stdout(args.output_file) as fout:
        fout.write(content)

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Count reads in BAM files')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('sequential_mapping', 
        help='generate Snakefile for sequential mapping')
    parser.add_argument('--rna-types', type=str, required=True,
        help='comma-separated list of rna types')
    parser.add_argument('--aligner', type=str, default='bowtie2',
        help='aligner to use')
    parser.add_argument('--output-file', '-o', type=str, default='-',
        help='output file')
    
    args = main_parser.parse_args()
    if args.command is None:
        main_parser.print_help()
        sys.exit(1)
    logger = logging.getLogger('generate_snakemake.' + args.command)

    command_handlers.get(args.command)(args)