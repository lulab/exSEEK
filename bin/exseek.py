#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
import yaml
import shutil
import shlex
import subprocess

def check_inputs(configfile):
    # check data_dir
    pass

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

    parser.add_argument('--step', type=str, required=True, choices=('quality_control', 
        'mapping', 'count_matrix', 'call_domains'))
    parser.add_argument('--dataset', '-d', type=str, required=True,
        help='dataset name')
    parser.add_argument('--config-dir', '-c', type=str, default='config',
        help='directory for configuration files')
    parser.add_argument('--cluster', action='store_true', help='submit to cluster')
    parser.add_argument('--cluster-config', type=str, 
        help='cluster configuration file ({config_dir}/cluster.yaml by default)')
    parser.add_argument('--cluster-command', type=str,
        help='command for submitting job to cluster (default read from {config_dir}/cluster_command.txt')
    parser.add_argument('--jobs', '-j', type=int, help='number of jobs')
    parser.add_argument('--dryrun', action='store_true', default=False, 
        help='do not actually run snakemake')
    parser.add_argument('--snakemake-args', type=str, 
        help='other snakemake arguments')
    args = parser.parse_args()

    logger = logging.getLogger('exseek')

    snakefile = None
    root_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    logger.info('root directory: {}'.format(root_dir))

    # check snakemake
    snakemake_path = shutil.which('snakemake')
    if snakemake_path is None:
        raise ValueError('cannot find snakemake command')
    
    # snakemake command
    snakemake_args = ['snakemake', '-k', '--rerun-incomplete']
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
    if args.step == 'quality_control':
        snakefile = os.path.join(root_dir, 'snakemake', 'quality_control.snakemake')
    elif args.step == 'mapping':
        if config['small_rna']:
            snakefile = os.path.join(root_dir, 'snakemake', 'mapping_small.snakemake')
        else:
            snakefile = os.path.join(root_dir, 'snakemake', 'mapping_long.snakemake')
    elif args.step == 'count_matrix':
        snakefile = os.path.join(root_dir, 'snakemake', 'count_matrix.snakemake')
    elif args.step == 'call_domains':
        snakefile = os.path.join(root_dir, 'snakemake', 'call_domains.snakemake')
    elif args.step == 'normalization':
        snakefile = os.path.join(root_dir, 'snakemake', 'normalization.snakemake')
    elif args.step == 'feature_selection':
        snakefile = os.path.join(root_dir, 'snakemake', 'feature_selection.snakemake')
    snakemake_args += ['--snakefile', snakefile, '--configfile', configfile]
    # number of jobs
    if args.jobs is not None:
        snakemake_args += ['-j', args.jobs]
    if args.dryrun is not None:
        snakemake_args += ['--dryrun']
    # extra args
    if args.snakemake_args is not None:
        snakemake_args += shlex.split(args.snakemake_args)
    # run snakemake
    snakemake_args = [str(s) for s in snakemake_args]
    logger.info('run snakemake: {}'.format(quoted_string_join(snakemake_args)))

    subprocess.check_call(snakemake_args, shell=False)
    