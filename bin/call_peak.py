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
def refine_peak_boundary(args):
    pass
    
if __name__ == '__main__':
    main_parser = argparse.ArgumentParser(description='Call peaks from exRNA signals')
    subparsers = main_parser.add_subparsers(dest='command')
    
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

    args = main_parser.parse_args()
    if args.command is None:
        raise ValueError('empty command')
    logger = logging.getLogger('call_peak.' + args.command)

    command_handlers.get(args.command)(args)

