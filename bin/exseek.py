#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
import yaml
import shutil
import shlex
import subprocess

def quoted_string_join(strs, sep=' '):
    quoted = []
    for s in strs:
        if len(s.split()) > 1:
            quoted.append('"' + s + '"')
        else:
            quoted.append(s)
    return sep.join(quoted)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='exSeek main program')

    parser.add_argument('step', type=str, 
        choices=('quality_control', 'quality_control_clean',
        'fastq_to_fasta', 'prepare_genome', 'bigwig',
        'mapping', 'count_matrix', 'call_domains', 'combine_domains',
        'normalization', 'feature_selection', 
        'differential_expression', 'evaluate_features', 'igv',
        'update_sequential_mapping', 'update_singularity_wrappers')
    )
    parser.add_argument('--dataset', '-d', type=str, required=True,
        help='dataset name')
    parser.add_argument('--config-dir', '-c', type=str, default='config',
        help='directory for configuration files')
    parser.add_argument('--cluster', action='store_true', help='submit to cluster')
    parser.add_argument('--cluster-config', type=str, 
        help='cluster configuration file ({config_dir}/cluster.yaml by default)')
    parser.add_argument('--cluster-command', type=str,
        help='command for submitting job to cluster (default read from {config_dir}/cluster_command.txt')
    parser.add_argument('--singularity', action='store_true',
        help='use singularity')
    args, extra_args = parser.parse_known_args()

    logger = logging.getLogger('exseek')

    snakefile = None
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    logger.info('root directory: {}'.format(root_dir))

    logger.info('read default config file')
    with open('snakemake/default_config.yaml', 'r') as f:
        default_config = yaml.load(f)

    # find snakemake executable
    snakemake_path = shutil.which('snakemake')
    if snakemake_path is None:
        raise ValueError('cannot find snakemake command')

    # snakemake command
    snakemake_args = ['snakemake', '-k', '--rerun-incomplete']
    extra_config = {}
    # check configuration file
    if not os.path.isdir(args.config_dir):
        raise ValueError('cannot find configuration directory: {}'.format(args.config_dir))
    configfile = os.path.join(args.config_dir, '{}.yaml'.format(args.dataset))
    if not os.path.isfile(configfile):
        raise ValueError('cannot find configuration file: {} '.format(configfile))
    with open(configfile, 'r') as f:
        config = yaml.load(f)
    # check cluster configuration
    if args.cluster:
        cluster_config = os.path.join(args.config_dir, 'cluster.yaml')
        if not os.path.isfile(cluster_config):
            if args.cluster_config is None:
                raise ValueError('cannot find {}/cluster.yaml and --cluster-config is not specified'.format(args.config_dir))
            cluster_config = args.cluster_config

        cluster_command_file = os.path.join(args.config_dir, 'cluster_command.txt')
        if os.path.isfile(cluster_command_file):
            with open(cluster_command_file, 'r') as f:
                cluster_command = f.read().strip()
        else:
            if args.cluster_command is None:
                raise ValueError('cannot find {}/cluster_command.txt and --cluster-command is not specified'.format(args.config_dir))
            cluster_command = args.cluster_command
        snakemake_args += ['--cluster', cluster_command, '--cluster-config', cluster_config]
    
    def update_sequential_mapping():
        snakefile = os.path.join(root_dir, 'snakemake', 'sequential_mapping.snakemake')
        logger.info('generate sequential_mapping.snakemake')
        subprocess.check_call([os.path.join(root_dir, 'bin/generate_snakemake.py'), 'sequential_mapping',
                '--rna-types', ','.join(config['rna_types']), 
                '-o', snakefile], shell=False)
        
    def update_singularity_wrappers():
        singularity_path = shutil.which('singularity')
        if singularity_path is None:
            raise ValueError('cannot find singularity')
        logger.info('generate singularity wrappers')
        subprocess.check_call([os.path.join(root_dir, 'bin/make_singularity_wrappers.py'), 
            '--image', default_config['singularity']['image'],
            '--list-file', os.path.join(root_dir, 'singularity/exports.txt'),
            '--singularity-path', singularity_path,
            '-o', default_config['singularity']['wrapper_dir']
        ], shell=False)
        
    # find proper version of snakemake
    if args.step == 'quality_control':
        if config['paired_end']:
            snakefile = os.path.join(root_dir, 'snakemake', 'quality_control_pe.snakemake')
        else:
            snakefile = os.path.join(root_dir, 'snakemake', 'quality_control_se.snakemake')
    elif args.step == 'quality_control_clean':
        if config['paired_end']:
            snakefile = os.path.join(root_dir, 'snakemake', 'quality_control_clean_pe.snakemake')
        else:
            snakefile = os.path.join(root_dir, 'snakemake', 'quality_control_clean_se.snakemake')
    elif args.step == 'mapping':
        if config['small_rna']:
            if not os.path.exists(os.path.join(root_dir, 'snakemake', 'sequential_mapping.snakemake')):
                update_sequential_mapping()
            snakefile = os.path.join(root_dir, 'snakemake', 'mapping_small.snakemake')
        else:
            if config['paired_end']:
                snakefile = os.path.join(root_dir, 'snakemake', 'mapping_long_pe.snakemake')
            else:
                snakefile = os.path.join(root_dir, 'snakemake', 'mapping_long_se.snakemake')
    elif args.step == 'count_matrix':
        if config['small_rna']:
            snakefile = os.path.join(root_dir, 'snakemake', 'count_matrix_small.snakemake')
        else:
            snakefile = os.path.join(root_dir, 'snakemake', 'count_matrix_long.snakemake')
    elif args.step == 'combine_domains':
        if config['small_rna']:
            snakefile = os.path.join(root_dir, 'snakemake', 'combine_domains_with_small.snakemake')
        else:
            raise ValueError('combine_domains can only be applied to small RNA-seq data')
    elif args.step == 'update_sequential_mapping':
        if config['small_rna']:
            update_sequential_mapping()
        sys.exit(0)
    elif args.step == 'update_singularity_wrappers':
        if args.singularity is None:
            raise ValueError('argument --singularity is required for step: update-singularity-wrappers')
        update_singularity_wrappers()
        sys.exit(0)
    elif args.step == 'bigwig':
        if config['small_rna']:
            snakefile = os.path.join(root_dir, 'snakemake', 'bigwig_small.snakemake')
        else:
            snakefile = os.path.join(root_dir, 'snakemake', 'bigwig_long.snakemake')
    elif args.step == 'call_domains':
        if config['small_rna']:
            snakefile = os.path.join(root_dir, 'snakemake', 'call_domains.snakemake')
        else:
            raise ValueError('call_domains can only be applied to small RNA-seq data')
    else:
        snakefile = os.path.join(root_dir, 'snakemake', args.step + '.snakemake')
    snakemake_args += ['--snakefile', snakefile, '--configfile', configfile]
    # set root_dir and bin_dir
    extra_config['bin_dir'] = os.path.join(root_dir, 'bin')
    extra_config['root_dir'] = root_dir
    extra_config['dataset'] = args.dataset
    # extra args
    snakemake_args = [str(s) for s in snakemake_args]
    snakemake_args += extra_args

    if args.singularity:
        if not os.path.isdir(default_config['singularity']['wrapper_dir']):
            update_singularity_wrappers()
        logger.info('enable singularity')
        extra_config['use_singularity'] = True
    
    # extra config
    snakemake_args += ['--config'] + ['{}={}'.format(key, val) for key, val in extra_config.items()]
    #subprocess.check_call(snakemake_args, shell=False)
    logger.info('run snakemake: {}'.format(quoted_string_join(snakemake_args)))
    # run snakemake
    os.execv(snakemake_path, snakemake_args)
    