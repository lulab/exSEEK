#! /usr/bin/env python

import sys
import argparse
import os
import pandas as pd

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

formats = ('table', 'csv', 'json', 'html', 'excel', 'hdf', 'sql', 'pickle')

def detect_format(filename):
    ext = os.path.splitext(filename)[1]
    if not ext:
        return
    format = {
        '.txt': 'table',
        '.csv': 'csv',
        '.json': 'json',
        '.xls': 'excel',
        '.xlsx': 'excel',
        '.h5': 'hdf',
        '.hdf5': 'hdf',
        '.hdf': 'hdf',
        '.sql': 'sql',
        '.pkl': 'pickle',
        '.pickle': 'pickle'
    }.get(ext)
    if format is None:
        raise ValueError('unknown file extension: {}'.format(ext))
    return format
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser('Convert table formats')
    parser.add_argument('--input-file', '-i', type=str, default='-')
    parser.add_argument('--output-file', '-o', type=str, default='-')
    parser.add_argument('--sformat', '-s', type=str,
        choices=formats, help='input format')
    parser.add_argument('--dformat', '-d', type=str, 
        choices=formats, help='output format')
    parser.add_argument('--reader-args', '-r', type=str, action='append')
    parser.add_argument('--writer-args', '-w', type=str, action='append')
    args = parser.parse_args()

    reader_args = {}
    if args.reader_args:
        for arg in args.reader_args:
            c = arg.split('=')
            if len(c) != 2:
                raise ValueError('reader args should be specified as key=value')

    sformat = detect_format(args.input_file)
    if not sformat:
        sformat = args.sformat
    if not sformat:
        raise ValueError('cannot detect format from input filename and --sformat is not specified')
    
    import pandas as pd
    read_df = {
        'table': pd.read_table,
        'csv': pd.read_csv,
        'json': pd.read_json,
        'html': pd.read_html,
        'excel': pd.read_excel,
        'hdf': pd.read_hdf,
        'sql': pd.read_sql,
        'pickle': pd.read_pickle
    }[sformat]
    
    with open_file_or_stdin(args.input_file) as f:
        df = read_df(f, **reader_args)

    dformat = detect_format(args.output_file)
    if not dformat:
        dformat = args.dformat
    if not dformat:
        raise ValueError('cannot detect format from output filename and --dformat is not specified')
    
    write_df = {
        'table': df.to_csv,
        'csv': df.to_csv,
        'json': df.to_json,
        'html': df.to_html,
        'excel': df.to_excel,
        'hdf': df.to_hdf,
        'sql': df.to_sql,
        'pickle': df.to_pickle
    }[dformat]

    writer_args = {}
    if args.writer_args:
        for arg in args.writer_args:
            c = arg.split('=')
            if len(c) != 2:
                raise ValueError('writer args should be specified as key=value')
            writer_args[c[0]] = c[1]

    with open_file_or_stdout(args.output_file) as f:
        write_df(f, **writer_args)
