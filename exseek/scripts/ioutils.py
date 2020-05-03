import os
import sys

def make_dir(path):
    if not os.path.exists(path):
        try:
            os.makedirs(path)
        except OSError as e:
            if e.errno == 17:
                pass
            else:
                raise e

def prepare_output_file(path):
    if path.startswith('/dev'):
        return
    dirpath = os.path.dirname(path)
    if dirpath == '':
        dirpath = '.'
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

def append_extra_line(f):
    """Yield an empty line after the last line in the file
    """
    for line in f:
        yield line
    yield ''

def open_file_or_stdout(filename):
    if filename == '-':
        return sys.stdout
    else:
        return open(filename, 'w')

def open_file_or_stdin(filename):
    if filename == '-':
        return sys.stdin
    else:
        return open(filename, 'r')

import zipfile
class ArchiveFile(object):
    def __init__(self, filename, mode='r', format='directory', **kwargs):
        self.filename = filename
        self.format = format
        if self.format == 'directory':
            if 'r' in mode:
                if not os.path.isdir(filename):
                    raise IOError('cannot open the directory: {}'.format(filename))
            else:
                os.makedirs(filename, exist_ok=True)
        elif self.format == 'file':
            self.f = open(file, mode)
        elif self.format == 'zip':
            if ('w' in mode) and ('compression' not in kwargs):
                kwargs['compression'] = zipfile.ZIP_STORED
            self.f = zipfile.ZipFile(filename, mode, **kwargs)

    def open(self, name, mode='r', **kwargs):
        if self.format == 'directory':
            return open(os.path.join(self.filename, name), mode)
        elif self.format == 'file':
            return self.f
        elif self.format == 'zip':
            if 'r' in mode:
                return self.f.open(name, mode)
            else:
                return

    def close(self):
        if self.format == 'directory':
            pass
        elif self.format == 'file':
            self.f.close()
        elif self.format == 'zipfile':
            self.f.close()

    def namelist(self):
        if self.format == 'directory':
            return os.listdir(self.filename)
        elif self.format == 'file':
            return self.filename
        elif self.format == 'zip':
            return self.f.namelist()
