from collections import namedtuple

GFFRecord = namedtuple('GFFRecord', ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attr'])

def read_gff(fin, zero_based=False):
    '''Read a GFF file
    Args:
        fin: file-like object
        zero_based: subtract start position by 1 if True
    Yields:
        feature: GFFRecord object
    '''
    for lineno, line in enumerate(fin):
        if line.startswith('#'):
            continue
        c = line.strip().split('\t')
        if len(c) < 9:
            raise ValueError('less than 9 columns found in GFF file at line{}'.format(lineno + 1))
        c[3] = int(c[3])
        c[4] = int(c[4])
        if zero_based:
            c[3] -= 1
        attrs = {}
        for a in c[8].split(';'):
            i = a.find('=')
            key = a[:i]
            val = a[(i + 1):]
            attrs[key] = val
        c[8] = attrs
        yield GFFRecord._make(c)