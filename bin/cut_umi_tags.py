#! /usr/bin/env python
import argparse, sys, os, errno
import re

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Remove UMI tags')
    parser.add_argument('--cut', '-u', type=int, default=0,
        help='trim n leading bases from 5\'-end')
    parser.add_argument('--input-file', '-i', type=str, default='-')
    parser.add_argument('--output-file', '-o', type=str, default='-')
    args = parser.parse_args()

    fin = sys.stdin
    if args.input_file != '-':
        fin = open(args.input_file, 'r')
    fout = sys.stdout
    if args.output_file != '-':
        fout = open(args.output_file, 'w')
    
    cut = args.cut
    record = []
    pat = re.compile('^[G]+')
    def process_record(record):
        record[1] = pat.sub('', record[1][cut:])
        record[3] = record[3][(len(record[3]) - len(record[1])):]
        fout.write('\n'.join(record))
        fout.write('\n')
    
    try:
        for i, line in enumerate(fin):
            if i % 4 == 0:
                if record:
                    process_record(record)
                    record = []
            record.append(line.strip())
        if record:
            process_record(record)
    except BrokenPipeError:
        pass
    except KeyboardInterrupt:
        pass
    finally:
        fin.close()
        fout.close()
