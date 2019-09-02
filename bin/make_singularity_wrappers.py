#! /usr/bin/env python
import json, os, sys
import argparse
import re
import subprocess
from io import StringIO
import stat
from glob import glob

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='get file list in conda packages')
    parser.add_argument('--image', type=str, required=True, 
        help='singularity image')
    parser.add_argument('--list-file', '-l', type=str, required=True,
        help='package list file')
    parser.add_argument('--output-dir', '-o', type=str, required=True,
        help='output directory')
    parser.add_argument('--singularity-path', type=str, default='singularity')
    parser.add_argument('--conda-path', type=str, default='/opt/conda')
    parser.add_argument('--runner', type=str, default='singularity',
        help='executable to run the container')
    parser.add_argument('--runner-path', type=str,
        help='path of the runner executable')
    args = parser.parse_args()

    is_inside_container = False
    for v in  ('SINGULARITY_CONTAINER', 'container_root'):
        if v in os.environ:
            is_inside_container = True
            break
    runner_path = args.runner
    if args.runner_path is not None:
        runner_path = args.runner_path
    # execute inside container
    if not is_inside_container:
        if args.runner == 'singularity':
            subprocess.check_call([args.runner, 'exec', args.image] + sys.argv, shell=False)
        elif args.runner == 'udocker':
            subprocess.check_call([args.runner, 'run', args.image] + sys.argv, shell=False)
        else:
            raise ValueError('unsupported runner')
        sys.exit(0)

    def make_wrapper(filename):
        #print('Make wrapper: {}'.format(filename))
        if not os.path.isdir(args.output_dir):
            os.makedirs(args.output_dir)
        basename = os.path.basename(filename)
        with open(os.path.join(args.output_dir, basename), 'w') as f:
            f.write('''#! /bin/bash
exec "{0}" exec "{1}" "{2}" "$@"
'''.format(runner_path, os.path.abspath(args.image), filename))
        os.chmod(os.path.join(args.output_dir, basename), 0o755)

    with open(args.list_file, 'r') as fin:
        for line in fin:
            c = line.strip().split()
            if len(c) != 2:
                continue
            pkg, source = c
            # find executable files in a apt package (Ubuntu or Debian OS)
            if source == 'apt':
                for filename in StringIO(unicode(subprocess.check_output(['dpkg', '-L', pkg], shell=False), encoding='utf-8')):
                    filename = filename.strip()
                    basename = os.path.basename(filename)
                    if os.path.basename(os.path.dirname(filename)) == 'bin':
                        st = os.stat(filename)
                        if (st.st_mode & stat.S_IXUSR) and stat.S_ISREG(st.st_mode):
                            make_wrapper(filename)
            # find executable files in a conda package
            elif source == 'conda':
                metadata_file = glob(os.path.join(args.conda_path, 'conda-meta', pkg + '-*.json'))
                if not metadata_file:
                    print('Warning: cannot find package: ' + pkg)
                    continue
                with open(metadata_file[0], 'r') as f:
                    metadata = json.load(f)
                    for filename in metadata['files']:
                        if os.path.dirname(filename) == 'bin':
                            make_wrapper(os.path.join(args.conda_path, filename))
    