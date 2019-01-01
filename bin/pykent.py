import ctypes
from ctypes.util import find_library
from ctypes import *
from enum import Enum
import sys, os
import numpy as np
import logging
logger = logging.getLogger('pykent')

__all__ = ['BigWigFile', 'init_library']

'''Python wrapper for UCSC kent library
=======================================

Steps to build the jkweb shared library:
* Download the source from http://hgdownload.soe.ucsc.edu/admin/exe/userApps.src.tgz
* Enter the kent/src directory.
* Make sure that openssl-devel, libpng-devel, zlib-devel, libuuid-devel are installed
* Find static library paths
for lib in png crypto ssl uuid;do
ld --verbose -Bstatic -l$lib 2>/dev/null | awk '{match($0, /attempt to open (.+) succeeded/,a); if(length(a)>0) print a[1]}'
done
* Modify lib/makefile to add the following lines:
%.o: %.c
        ${CC} -fPIC ${COPT} ${CFLAGS} ${HG_DEFS} ${LOWELAB_DEFS} ${HG_WARN} ${HG_INC} ${XINC} -o $@ -c $<

$(MACHTYPE)/libjkweb.so: $(O) $(MACHTYPE)
    gcc -fPIC -shared -o $(MACHTYPE)/libjkweb.so $(O) -Wl,-z,defs -L../htslib -lhts -lm -lz -lpthread -lpng -lcrypto -lssl -luuid
* Enter the htslib directory, run the following commands to generate libhts.a
CFLAGS="-fPIC -DUCSC_CRAM=0 -DKNETFILE_HOOKS=1" ./configure
make
* Enter the lib/ directory. Run
make x86_64/libjkweb.so
'''


'''Definitions
======================
inc/common.h
#define UBYTE unsigned char   /* Wants to be unsigned 8 bits. */
#define BYTE signed char      /* Wants to be signed 8 bits. */
#define UWORD unsigned short  /* Wants to be unsigned 16 bits. */
#define WORD short	      /* Wants to be signed 16 bits. */
#define bits64 unsigned long long  /* Wants to be unsigned 64 bits. */
#define bits32 unsigned       /* Wants to be unsigned 32 bits. */
#define bits16 unsigned short /* Wants to be unsigned 16 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */
#define signed32 int	      /* Wants to be signed 32 bits. */
#define bits8 unsigned char   /* Wants to be unsigned 8 bits. */
inc/bits.h
typedef unsigned char Bits;
'''

UBYTE = c_ubyte
BYTE = c_byte
UWORD = c_ushort
WORD = c_short
bits64 = c_ulonglong
bits32 = c_uint
bits16 = c_ushort
bits8 = c_ubyte
boolean = c_int
Bits = c_ubyte

'''Bits
#define bitToByteSize(bitSize) ((bitSize+7)/8)
'''

def bitToByteSize(bitSize):
    return (bitSize+7)//8

'''
struct bbiFile 
/* An open bbiFile */
    {
    struct bbiFile *next;	/* Next in list. */
    char *fileName;		/* Name of file - for better error reporting. */
    struct udcFile *udc;	/* Open UDC file handle. */
    bits32 typeSig;		/* bigBedSig or bigWigSig for now. */
    boolean isSwapped;		/* If TRUE need to byte swap everything. */
    struct bptFile *chromBpt;	/* Index of chromosomes. */
    bits16 version;		/* Version number - initially 1. */
    bits16 zoomLevels;		/* Number of zoom levels. */
    bits64 chromTreeOffset;	/* Offset to chromosome index. */
    bits64 unzoomedDataOffset;	/* Start of unzoomed data. */
    bits64 unzoomedIndexOffset;	/* Start of unzoomed index. */
    bits16 fieldCount;		/* Number of columns in bed version. */
    bits16 definedFieldCount;   /* Number of columns using bed standard definitions. */
    bits64 asOffset;		/* Offset to embedded null-terminated AutoSQL file. */
    bits64 totalSummaryOffset;	/* Offset to total summary information if any.  
				   (On older files have to calculate) */
    bits32 uncompressBufSize;	/* Size of uncompression buffer, 0 if uncompressed */
    bits64 extensionOffset;	/* Start of header extension block or 0 if none. */
    struct cirTreeFile *unzoomedCir;	/* Unzoomed data index in memory - may be NULL. */
    struct bbiZoomLevel *levelList;	/* List of zoom levels. */

    /* Fields based on extension block. */
    bits16 extensionSize;   /* Size of extension block */
    bits16 extraIndexCount; /* Number of extra indexes (on fields other than chrom,start,end */ 
    bits64 extraIndexListOffset;    /* Offset to list of extra indexes */
    };
'''

class udcFile(Structure):
    pass

class bptFile(Structure):
    pass

class circTreeFile(Structure):
    pass

'''struct bbiZoomLevel
/* A zoom level in bigWig file. */
    {
    struct bbiZoomLevel *next;		/* Next in list. */
    bits32 reductionLevel;		/* How many bases per item */
    bits32 reserved;			/* Zero for now. */
    bits64 dataOffset;			/* Offset of data for this level in file. */
    bits64 indexOffset;			/* Offset of index for this level in file. */
    };
'''
class bbiZoomLevel(Structure):
    pass

bbiZoomLevel._fields_ = [
    ('next', POINTER(bbiZoomLevel)),
    ('reductionLevel', bits32),
    ('reserved', bits32),
    ('dataOffset', bits64),
    ('indexOffset', bits64)
]

class bbiFile(Structure):
    pass

class bbiChromIdSize(Structure):
    _fields_ = [
        ('chromId', c_uint),
        ('chromSize', c_uint)
    ]

class bbiChromInfo(Structure):
    pass

bbiChromInfo._fields_ = [
    ('next', POINTER(bbiChromInfo)),
    ('name', c_char_p),
    ('id', c_uint),
    ('size', c_uint)
]

class bbiInterval(Structure):
    pass

bbiInterval._fields_ = [
    ('next', POINTER(bbiInterval)),
    ('start', c_uint),
    ('end', c_uint),
    ('val', c_double)
]

'''
enum bbiSummaryType
/* Way to summarize data. */
    {
    bbiSumMean = 0,	/* Average value */
    bbiSumMax = 1,	/* Maximum value */
    bbiSumMin = 2,	/* Minimum value */
    bbiSumCoverage = 3,  /* Bases in region containing actual data. */
    bbiSumStandardDeviation = 4, /* Standard deviation in window. */
    };
'''
bbiSumMean = 0
bbiSumMax = 1
bbiSumMin = 2
bbiSumCoverage = 3
bbiSumStandardDeviation = 4

'''
struct bbiSummary
/* A summary type item. */
    {
    struct bbiSummary *next;
    bits32 chromId;		/* ID of associated chromosome. */
    bits32 start,end;		/* Range of chromosome covered. */
    bits32 validCount;		/* Count of (bases) with actual data. */
    float minVal;		/* Minimum value of items */
    float maxVal;		/* Maximum value of items */
    float sumData;		/* sum of values for each base. */
    float sumSquares;		/* sum of squares for each base. */
    bits64 fileOffset;		/* Offset of summary in file. */
    };
'''
class bbiSummary(Structure):
    pass

bbiSummary._fields_ = [
    ('next', POINTER(bbiSummary)),
    ('chromId', bits32),
    ('start', bits32),
    ('end', bits32),
    ('validCount', bits32),
    ('minVal', c_float),
    ('maxVal', c_float),
    ('sumData', c_float),
    ('sumSquares', c_float),
    ('fileOffset', bits64)
]

'''
struct bbiSummaryElement
/* An element of a summary from the user side. */
    {
    bits64 validCount;		/* Count of (bases) with actual data. */
    double minVal;		/* Minimum value of items */
    double maxVal;		/* Maximum value of items */
    double sumData;		/* sum of values for each base. */
    double sumSquares;		/* sum of squares for each base. */
    };
'''

class bbiSummaryElement(Structure):
    _fields_ = [
        ('validCount', bits64),
        ('minVal', c_double),
        ('maxVal', c_double),
        ('sumData', c_double),
        ('sumSquares', c_double)
    ]

class lm(Structure):
    pass

'''
struct bigWigValsOnChrom
/* Object for bulk access a chromosome at a time.  This is faster than
 * doing bigWigInterval queries when you have ~3000 or more queries. */
     {
     struct bigWigValsOnChrom *next;
     char *chrom;	/* Current chromosome. */
     long chromSize;	/* Size of current chromosome. */
     long bufSize;	/* Size of allocated buffer */
     double *valBuf;	/* A value for each base on chrom. Zero where no data. */
     Bits *covBuf;	/* A bit for each base with data. */
     };
'''
class bigWigValsOnChrom(Structure):
    pass

bigWigValsOnChrom._fields_ = [
    ('next', POINTER(bigWigValsOnChrom)),
    ('chrom', c_char_p),
    ('chromSize', c_long),
    ('bufSize', c_long),
    ('valBuf', POINTER(c_double)),
    ('covBuf', POINTER(Bits))
]

'''Standard C library
#include <stdio.h>

extern FILE *stdin;
extern FILE *stdout;
extern FILE *stderr;
'''
class FILE(Structure):
    pass

FILE_p = POINTER(FILE)


class SharedLibrary(object):
    def __init__(self, name):
        if os.path.isfile(name):
            libpath = name
        else:
            libpath = find_library(name)
            if libpath is None:
                raise ValueError('cannot find library {}'.format(name))
        self.lib = cdll.LoadLibrary(libpath)
        self.functions = {}
        self.variables = {}
    
    def resolve_function(self, name, restype, argtypes):
        '''
        '''
        try:
            f = getattr(lib, name)
            f.restype = restype
            if not isinstance(argtypes, tuple):
                if getattr(argtypes, '__iter__'):
                    argtypes = tuple(argtypes)
                else:
                    argtypes = (argtypes,)
            f.argtypes = argtypes
            self.functions[name] = f
        except AttributeError:
            logger.warn('failed to resolve function {}'.format(name))
    
    def resolve_variable(self, name, vartype):
        try:
            v = vartype.in_dll(self.lib, name)
            self.variables[name] = v
        except AttributeError:
            logger.warn('failed to resolve variable {}'.format(name))
    
    def export_to_globals(self):
        '''Export all resolved symbols to global variables
        '''
        for name, value in self.functions.items():
            globals()[name] = value
        for name, value in self.variables.items():
            globals()[name] = value

def init_library(libpath='jkweb'):
    '''Load shared library and resolve symbols
    libpath: path to libjkweb.so. Find in LD_LIBRARY_PATH if libpath is None.
    '''
    lib = SharedLibrary(libpath)

    # struct lm *lmInit(int blockSize);
    lib.resolve_function('lmInit', POINTER(lm), (c_int,))
    # void lmCleanup(struct lm **pLm);
    lib.resolve_function('lmCleanup', None, (POINTER(POINTER(lm)),))
    # struct bbiFile *bigWigFileOpen(char *fileName);
    lib.resolve_function('bigWigFileOpen', POINTER(bbiFile), (c_char_p,))
    # struct bbiFile *bbiFileOpen(char *fileName, bits32 sig, char *typeName);
    lib.resolve_function('bbiFileOpen', POINTER(bbiFile), (c_char_p, bits32, c_char_p))
    # void bbiFileClose(struct bbiFile **pBwf);
    lib.resolve_function('bbiFileClose', None, (POINTER(POINTER(bbiFile)),))
    # define bigWigFileClose(a) bbiFileClose(a)
    bigWigFileClose = bbiFileClose
    # bits32 bbiChromSize(struct bbiFile *bbi, char *chrom);
    lib.resolve_function('bbiChromSize', bits32, (POINTER(bbiFile), c_char_p))
    # struct bbiChromInfo *bbiChromList(struct bbiFile *bbi);
    lib.resolve_function('bbiChromList', POINTER(bbiChromInfo), (POINTER(bbiFile),))
    # void bbiChromInfoFreeList(struct bbiChromInfo **pList);
    lib.resolve_function('bbiChromInfoFreeList', None, (POINTER(POINTER(bbiChromInfo)),))
    # struct bbiInterval *bigWigIntervalQuery(struct bbiFile *bwf, char *chrom, bits32 start, bits32 end,
	# struct lm *lm);
    lib.resolve_function('bigWigIntervalQuery', POINTER(bbiInterval), 
        (POINTER(bbiFile), c_char_p, bits32, bits32, POINTER(lm)))
    # boolean isBigWig(char *fileName);
    lib.resolve_function('isBigWig', boolean, (c_char_p,))
    # boolean bigWigFileCheckSigs(char *fileName);
    lib.resolve_function('bigWigFileCheckSigs', boolean, (c_char_p,))
    # struct bigWigValsOnChrom *bigWigValsOnChromNew();
    lib.resolve_function('bigWigValsOnChromNew', POINTER(bigWigValsOnChrom), ())
    # void bigWigValsOnChromFree(struct bigWigValsOnChrom **pChromVals);
    lib.resolve_function('bigWigValsOnChromFree', None, (POINTER(POINTER(bigWigValsOnChrom)),))
    # boolean bigWigValsOnChromFetchData(struct bigWigValsOnChrom *chromVals, char *chrom, 
	# struct bbiFile *bigWig);
    lib.resolve_function('bigWigValsOnChromFetchData', boolean, 
        (POINTER(bigWigValsOnChrom), c_char_p, POINTER(bbiFile)))
    # boolean bigWigSummaryArray(struct bbiFile *bwf, char *chrom, bits32 start, bits32 end,
	# enum bbiSummaryType summaryType, int summarySize, double *summaryValues);
    lib.resolve_function('bigWigSummaryArray', boolean, 
        (POINTER(bbiFile), c_char_p, bits32, bits32, c_int, c_int, POINTER(c_double)))
    # boolean bigWigSummaryArrayExtended(struct bbiFile *bwf, char *chrom, bits32 start, bits32 end,
	# int summarySize, struct bbiSummaryElement *summary);
    lib.resolve_function('bigWigSummaryArrayExtended', boolean,
        (POINTER(bbiFile), c_char_p, bits32, bits32, c_int, POINTER(bbiSummaryElement)))
    # double bigWigSingleSummary(struct bbiFile *bwf, char *chrom, int start, int end,
    # enum bbiSummaryType summaryType, double defaultVal);
    lib.resolve_function('bigWigSingleSummary', c_double,
        (POINTER(bbiFile), c_char_p, c_int, c_int, c_int, c_double))
    # boolean bitReadOne(Bits *b, int bitIx);
    lib.resolve_function('bitReadOne', boolean, (POINTER(Bits), c_int))
    
    # export symbols
    lib.export_to_globals()

    # resolve the C library
    lib = SharedLibrary('c')
    lib.resolve_variable('stdin', FILE_p)
    lib.resolve_variable('stdout', FILE_p)
    lib.resolve_variable('stderr', FILE_p)
    lib.resolve_function('setbuf', None, (FILE_p, c_char_p))
    lib.resolve_function('setbuffer', None, (FILE_p, c_char_p, c_size_t))
    

class BigWigFile(object):
    def __init__(self, filename):
        filename = bytes(filename, encoding='utf-8')
        self.bwf = bigWigFileOpen(filename)
        if not self.bwf:
            raise IOError('cannot open the bigwig file: ' + filename)
        self.lm = lmInit(0)

    def __del__(self):
        if self.bwf is not None:
            self.close()

    @staticmethod
    def is_bigwig(filename):
        filename = bytes(filename, encoding='utf-8')
        return bool(isBigWig(filename))
    
    @staticmethod
    def check_sigs(filename):
        filename = bytes(filename, encoding='utf-8')
        return bool(bigWigFileCheckSigs(filename))

    def close(self):
        bbiFileClose(self.bwf)
        lmCleanup(byref(self.lm))
        self.lm = None
        self.bwf = None
    
    def get_chrom_list(self):
        chrom_list = []
        chrom_info_head = bbiChromList(self.bwf)
        chrom_info = chrom_info_head
        while chrom_info:
            chrom_info = chrom_info.contents
            chrom_list.append((str(chrom_info.name, encoding='utf-8'), chrom_info.size))
            chrom_info = chrom_info.next
        bbiChromInfoFreeList(byref(chrom_info))
        return chrom_list

    def interval_query(self, chrom, start, end, dtype='float', fillna=None):
        """Returns a numpy array of all values in range chrom:start-end
        """
        chrom = bytes(chrom, encoding='utf-8')
        interval = bigWigIntervalQuery(self.bwf, chrom, start, end, self.lm)
        if fillna is None:
            fillna = np.nan
        values = np.full(end - start, fillna, dtype=dtype)
        while interval:
            interval = interval.contents
            values[(interval.start - start) : (interval.end - end)] = interval.val
            interval = interval.next
        return values
    
    def values_on_chrom(self, chrom, dtype='float', fillna=None):
        chrom = bytes(chrom, encoding='utf-8')
        p_chrom_vals = bigWigValsOnChromNew()
        if bigWigValsOnChromFetchData(p_chrom_vals, chrom, self.bwf):
            chrom_vals = p_chrom_vals.contents
            values = np.ctypeslib.as_array(chrom_vals.valBuf, shape=(chrom_vals.bufSize,)).astype(dtype)
            cov = np.ctypeslib.as_array(chrom_vals.covBuf, shape=(bitToByteSize(chrom_vals.bufSize),))
            #true_bits = np.asarray(unpackbits(chrom_vals.covBuf, chrom_vals.bufSize))
            cov = np.unpackbits(cov)[:chrom_vals.bufSize]
            #if not np.all(true_bits == cov):
            #    raise ValueError('incorrect unpackbits')
            if fillna is None:
                fillna = np.nan
            values[cov == 0] = fillna
        bigWigValsOnChromFree(byref(p_chrom_vals))
        return values

    def interval_query_blocked(self, chrom, start, end, block_starts, block_sizes):
        """Similar to interval_query except that only values in the blocks are fetched
        """
        n_blocks = len(block_starts)
        value_starts = np.cumsum([0] + block_sizes)[:-1]
        length = np.sum(block_sizes)
        values = np.full(length, np.nan, dtype='float64')

        for i, block_start, block_size in zip(range(n_blocks), block_starts, block_sizes):
            block_start += start
            interval = bigWigIntervalQuery(self.bwf, chrom,
                block_start, block_start + block_size, self.lm)
            while interval:
                interval = interval.contents
                values[(interval.start - block_start + value_starts[i]) :
                    (interval.end - block_start + value_starts[i])] = interval.val
                interval = interval.next
        return values

    def __del__(self):
        if self.bwf is not None:
            self.close()

def unpackbits(bits, bitCount):
    return [bitReadOne(bits, i) for i in range(bitCount)]

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser('BigWig library')
    parser.add_argument('--input-file', '-i', type=str, required=True, help='input BigWig file')
    parser.add_argument('--library', '-l', type=str, required=True, help='path to libjkweb.so')
    args = parser.parse_args()

    init_library(args.library)

    print('is_bigwig:', BigWigFile.is_bigwig(args.input_file))
    print('check_sigs:', BigWigFile.check_sigs(args.input_file))
    bwf = BigWigFile(args.input_file)
    chrom_sizes = bwf.get_chrom_list()
    print('chrom_list', chrom_sizes)
    #for chrom, size in chrom_sizes:
    #    values = bwf.interval_query(chrom, 0, size, na_rep=0).astype('int')
    #    print('interval_query {}:0-{}: {}'.format(chrom, size, ','.join(values.astype('str'))))
    for chrom, size in chrom_sizes:
        values = bwf.values_on_chrom(chrom, dtype='int', na_rep=0)
        print('interval_query {}:0-{}: {}'.format(chrom, size, ','.join(values.astype('str'))))