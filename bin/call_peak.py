#! /usr/bin/env python
import argparse, sys, os, errno
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')

import numpy as np
'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
sns.set()
'''
import pandas as pd
from pandas import DataFrame, Series
from scipy.fftpack import fft
from scipy.signal import convolve
import numba

command_handlers = {}
def command_handler(f):
    command_handlers[f.__name__] = f
    return f

def read_coverage(filename):
    coverage = []
    gene_ids = []
    with open(filename, 'r') as f:
        for line in f:
            c = line.strip().split('\t')
            gene_id = c[0]
            values = np.array(c[1:]).astype(np.float64)
            gene_ids.append(gene_id)
            coverage.append(values)
    return gene_ids, coverage

def plot_peak_calling(sig, gene_id, 
        local_bg_weight=0.5, bg_global=None,
        min_width=1, max_width=7):
    df = []
    positions = np.arange(len(sig))
    df.append(DataFrame({'Position': positions, 
                         'Signal': sig, 
                         'Name': 'Original (gene_id = {})'.format(gene_id),
                         'Color': np.full(len(sig), 0, dtype=np.int32)}))
    
    if bg_global is None:
        bg_global = np.mean(sig)

    for width in range(min_width, max_width + 1):
        filter = np.full(width, 1.0/width)
        bg_local = convolve(sig, filter, mode='same')
        bg = local_bg_weight*bg_local + (1.0 - local_bg_weight)*bg_global
        bg[np.isclose(bg, 0)] = 1
        snr = sig/bg
        peaks = (snr > 1.0).astype(np.int32)
        smoothed_peaks = icm_smooth(peaks, h=-2.0, beta=4.0, eta=2.0)
        df.append(DataFrame({'Position': positions,
                            'Signal': bg,
                            'Name': np.full(len(sig), 'background (width={})'.format(width)),
                            'Color': np.full(len(sig), 1, dtype=np.int32)}))
        df.append(DataFrame({'Position': positions,
                            'Signal': peaks,
                            'Name': np.full(len(sig), 'peak (width={})'.format(width)),
                            'Color': np.full(len(sig), 2, dtype=np.int32)}))
        df.append(DataFrame({'Position': positions,
                            'Signal': smoothed_peaks,
                            'Name': np.full(len(sig), 'smoothed peak (width={})'.format(width)),
                            'Color': np.full(len(sig), 3, dtype=np.int32)}))

    df = pd.concat(df, axis=0)
    g = sns.FacetGrid(df, row='Name', hue='Color', row_order=df['Name'].unique(), size=0.7, aspect=32, sharey=False, 
                    xlim=(0, max(200, len(sig))))
    #g.map(plt.plot, 'Position', 'Signal', drawstyle='steps', linewidth=1.5, fillstyle='full')
    g.map(plt.fill_between, 'Position', 'Signal', step='pre', edgecolor='none')
    g.fig.subplots_adjust(hspace=0.8)
    return g

@command_handler
def plot(args):
    from tqdm import tqdm

    logger.info('read input file: ' + args.input_file)
    gene_ids, signals = read_coverage(args.input_file)
    if args.use_log:
        signals = [np.log10(np.maximum(1e-3, a)) + 3 for a in signals]
    #logsignal_mean = np.array([np.mean(np.log10(np.maximum(a, 0.001))) for a in signals['exRNA_rpkm']])
    #logsignal_mean_nonzero = logsignal_mean[~np.isclose(logrpkm_mean, -3)]
    signals_mean = np.asarray([np.mean(a) for a in signals])
    bg_global = np.median(signals_mean)

    logger.info('create output plot file: ' + args.output_file)
    sns.set_style("dark", {'font.family': 'Arial'})
    with PdfPages(args.output_file) as pdf:
        if args.n_genes is not None:
            indices = np.random.choice(len(signals), size=args.n_genes, replace=False)
            indices = np.arange(args.n_genes)
        elif args.genes is not None:
            genes = args.genes.split(',')
            gene_id_dict = {gene_id:i for i, gene_id in enumerate(gene_ids)}
            indices = []
            for gene_id in genes:
                indices.append(gene_id_dict[gene_id])
            indices = np.array(indices)
        for i in tqdm(indices, unit='gene'):
            g = plot_peak_calling(signals[i], gene_ids[i], bg_global=bg_global, 
                min_width=args.min_width, max_width=args.max_width)
            pdf.savefig(g.fig)
            plt.close()

@numba.jit('int64(int32[:], int32[:], float64, float64, float64)')
def icm_update(x, y, h=0.0, beta=1.0, eta=2.1):
    n_changes = 0
    N = x.shape[0]
    for i in range(N):
        dx = -2*x[i]
        dE = 0
        if i > 0:
            dE += h*dx - beta*dx*x[i - 1] - eta*dx*y[i]
        if i < (N - 1):
            dE += h*dx - beta*dx*x[i + 1] - eta*dx*y[i]
        if dE < 0:
            x[i] = -x[i]
            n_changes += 1
    return n_changes
        
def icm_smooth(x, h=0.0, beta=1.0, eta=2.1):
    '''Smooth signals using iterated conditional modes
    Args:
        x: 1D signal
    Returns:
        Smoothed signal of the same length of x
    '''
    x = x*2 - 1
    y = x.copy()
    #E = h*np.sum(x) - beta*x[:-1]*x[1:] - eta*x*y
    n_updates = icm_update(x, y, h=h, beta=beta, eta=eta)
    while n_updates > 0:
        n_updates = icm_update(x, y, h=h, beta=beta, eta=eta)
    x = (x > 0).astype(np.int32)
    return x

def call_peak_gene(sig, local_bg_width=3, local_bg_weight=0.5, bg_global=None, smooth=False):
    if bg_global is None:
        bg_global = np.mean(sig)

    filter = np.full(local_bg_width, 1.0/local_bg_width)
    bg_local = convolve(sig, filter, mode='same')
    bg = local_bg_weight*bg_local + (1.0 - local_bg_weight)*bg_global
    bg[np.isclose(bg, 0)] = 1
    snr = sig/bg
    peaks = (snr > 1.0).astype(np.int32)
    if smooth:
        peaks = icm_smooth(peaks, h=-2.0, beta=4.0, eta=2.0)
    x = np.zeros(len(peaks) + 2, dtype=np.int32)
    x[1:-1] = peaks
    starts = np.nonzero(x[1:] > x[:-1])[0]
    ends = np.nonzero(x[:-1] > x[1:])[0]
    peaks = np.column_stack([starts, ends])
    return peaks

def estimate_bg_global(signals):
    '''signals
    '''
    signals_mean = np.asarray([np.mean(s) for s in signals])
    bg = np.median(signals_mean)
    return bg

def call_peaks(signals, min_length=None):
    bg_global = estimate_bg_global(signals)
    peaks = []
    for i, signal in enumerate(signals):
        peak_locations = call_peak_gene(signal, bg_global=bg_global, smooth=True)
        for start, end in peak_locations:
            if (min_length is not None) and ((end - start) >= min_length):
                peaks.append((i, start, end))
    return peaks

@command_handler
def call_peak(args):
    from tqdm import trange

    logger.info('read input file: ' + args.input_file)
    gene_ids, signals = read_coverage(args.input_file)
    if args.use_log:
        signals = [np.log10(np.maximum(1e-3, a)) + 3 for a in signals]

    signals_mean = np.asarray([np.mean(a) for a in signals])
    bg_global = np.median(signals_mean)

    logger.info('create output plot file: ' + args.output_file)
    with open(args.output_file, 'w') as fout:
        for i in trange(len(signals), unit='gene'):
            peaks_locations = call_peak_gene(signals[i], bg_global=bg_global, local_bg_weight=args.local_bg_weight,
                local_bg_width=args.local_bg_width, smooth=args.smooth)
            for j in range(peaks_locations.shape[0]):
                fout.write('{}\t{}\t{}\n'.format(gene_ids[i], peaks_locations[j, 0], peaks_locations[j, 1]))

@command_handler
def gviz(args):
    import string
    import subprocess

    r_template = string.Template('''library(Gviz)
options(ucscChromosomeNames=FALSE)
message('read peak file: ', '$peaks_file')
peakTrack <- AnnotationTrack('$peaks_file', 
                             genome='hg38', name='peak', shape='box')
message('read genomecov file: ', '$genomecov_plus_file')
covPlusTrack <- DataTrack('$genomecov_plus_file', genome='hg38', type='histogram',
    lty='blank', name='Read coverage (+)')
message('read genomecov file: ', '$genomecov_minus_file')
covMinusTrack <- DataTrack('$genomecov_minus_file', genome='hg38', type='histogram', 
    lty='blank', name='Read coverage (-)')
library(biomaRt)
bm <- useMart(host = "useast.ensembl.org",
              biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
grTrack <- BiomartGeneRegionTrack(biomart=bm, showFeatureId=TRUE, showId=TRUE, stacking='squish', genome='hg38', name='Genes')
gtrack <- GenomeAxisTrack()
plot_regions <- read.table('$region_file', header=FALSE)
pdf('$output_file', width=10, height=3)
for(i in 1:nrow(plot_regions)){
    message(sprintf('plot region %s (%s:%d-%d)', plot_regions[i, 4], plot_regions[i, 1], plot_regions[i, 2], plot_regions[i, 3]))
    plotTracks(list(grTrack, gtrack, covPlusTrack, covMinusTrack, peakTrack), sizes=c(3, 1, 1, 1, 1),
            chromosome=plot_regions[i, 1], from=plot_regions[i, 2], to=plot_regions[i, 3])
}
dev.off()
''')
    r_source = r_template.safe_substitute(peaks_file=args.peaks_file,
        genomecov_plus_file=args.genomecov_plus_file,
        genomecov_minus_file=args.genomecov_minus_file,
        region_file=args.region_file,
        output_file=args.output_file)
    p = subprocess.Popen(['R', '--vanilla', '--slave', '--no-save', '--no-restore'], stdin=subprocess.PIPE)
    p.communicate(bytes(r_source, encoding='ascii'))

if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Call peaks from exRNA signals')
    subparsers = main_parser.add_subparsers(dest='command')

    parser = subparsers.add_parser('plot')
    parser.add_argument('--input-file', '-i', type=str, required=True,
        help='input file of exRNA signals for each transcript')
    g = parser.add_mutually_exclusive_group(required=True)
    g.add_argument('--n-genes', type=int,
        help='number of genes to plot')
    g.add_argument('--genes', type=str, default='',
        help='comma separated list of gene ids')
    parser.add_argument('--use-log', action='store_true', 
        help='use log10 instead raw signals')
    parser.add_argument('--min-width', type=int, default=1)
    parser.add_argument('--max-width', type=int, default=7)
    parser.add_argument('--local-bg-weight', type=float, default=0.5, 
        help='weight for local background (0.0-1.0)')
    parser.add_argument('--output-file', '-o', type=str, required=True,
        help='output plot file (PDF format)')
    
    parser = subparsers.add_parser('call_peak')
    parser.add_argument('--input-file', '-i', type=str, required=True,
        help='input file of exRNA signals for each transcript')
    parser.add_argument('--use-log', action='store_true', 
        help='use log10 instead raw signals')
    parser.add_argument('--smooth', action='store_true',
        help='merge adjacent peaks')
    parser.add_argument('--local-bg-width', type=int, default=3,
        help='number of nearby bins for estimation of local background')
    parser.add_argument('--local-bg-weight', type=float, default=0.5, 
        help='weight for local background (0.0-1.0)')
    parser.add_argument('--output-file', '-o', type=str, required=True,
        help='output plot file BED format')

    parser = subparsers.add_parser('gviz')
    parser.add_argument('--peaks-file', type=str, required=True,
        help='peaks in BED format')
    parser.add_argument('--genomecov-plus-file', type=str, required=True,
        help='genome coverage in BigWig format (+ strand)')
    parser.add_argument('--genomecov-minus-file', type=str, required=True,
        help='genome coverage in BigWig format (- strand)')
    parser.add_argument('--region-file', type=str, required=True,
        help='specify region to plot in BED format')
    parser.add_argument('--output-file', '-o', type=str, required=True,
        help='output file in PDF format')
    
    parser = subparsers.add_parser('heatmap')
    parser.add_argument('--')

    args = main_parser.parse_args()
    if args.command is None:
        raise ValueError('empty command')
    logger = logging.getLogger('call_peak.' + args.command)

    command_handlers.get(args.command)(args)

