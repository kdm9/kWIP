# Copyright 2016 Kevin Murray <spam@kdmurray.id.au>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from __future__ import print_function, division, absolute_import

import numpy as np

import itertools as itl
from glob import glob
import os
from os import path
from sys import stderr, stdout
from multiprocessing import Pool
from functools import partial
import argparse
from argparse import (
    ArgumentParser,
    RawDescriptionHelpFormatter as RawFormatter,
)
from textwrap import dedent

from kmkm import KmerCollection, KmerCounter
from .logging import *


def countone(reads, ksize, cvsize, tables):
    kc = KmerCounter(ksize=ksize, cvsize=cvsize, cbf_tables=tables)
    kc.count_file(reads)
    return reads, kc.counts().copy()

def count_main():
    desc = "Counts k-mers into an aggregated kmer count matrix."
    epilog = dedent("""\
    Will use about 6 * CVLEN bytes of RAM per file (or 2x with --no-cms).
    """)
    #An optional pre-counting command for e.g. QC or SRA dumping can be given
    #with `-p`. The pre-command can be a shell pipeline combining the effects of
    #multiple programs. Interleaved or single ended reads must be printed on
    #stdout by the command(s). The pre-command uses the find/xargs/GNU Parallel
    #convention of using a pair of '{}' to mark where the filename should be
    #placed. Examples of a pre-command include:

    #    --precmd 'fastq-dump --split-spot --stdout {}'
    #    --precmd 'zcat {} | trimit'


    parser = ArgumentParser(description=desc, epilog=epilog,
                            formatter_class=RawFormatter)
    parser.add_argument(
        '-k', '--ksize', type=int, default=21,
        help='K-mer length')
    parser.add_argument(
        '-v', '--cvsize', type=float, default=5e8,
        help='Count vector length')
    parser.add_argument(
        '-j', '--jobs', default=1, type=int,
        help='Number of parallel jobs')
    #parser.add_argument(
    #    '-p', '--precmd', required=False,
    #    help='Shell pipeline to run on input files before hashing')
    parser.add_argument(
        '--no-cms', action='store_false', dest='use_cms',
        help='Disable the CMS counter, use only a count vector')
    parser.add_argument('outfile', help='Output file/directory')
    parser.add_argument('readfiles', nargs='+', help='Read files')
    args = parser.parse_args()

    info("Counting sequence files...")
    #mode = "a" if args.append else "w"
    mode = "w"
    tables = 3 if args.use_cms else 0
    counttable = KmerCollection(
        args.outfile, mode=mode, ksize=args.ksize, cvsize=args.cvsize,
        cbf_tables=tables)

    if args.jobs > 1:
        pool = Pool(args.jobs)
        mapper = pool.imap_unordered
    else:
        mapper = map

    f = partial(countone, ksize=args.ksize, cvsize=args.cvsize, tables=tables)
    for readfile, kc in mapper(f, args.readfiles):
        counttable.add_counts(kc, readfile)
        progress(readfile, kc.sum())
    info("All done!")


def weight_main():
    outfile = opts['WEIGHTFILE']
    countfiles = opts['COUNTFILES']

    popfreq = None
    nsamples = len(countfiles)

    for countfile in countfiles:
        info("Loading",  countfile, end='... ')
        counts = read_array(countfile, name='counts')
        if popfreq is None:
            popfreq = np.zeros(counts.shape, dtype=np.float32)
        popfreq_add_sample(popfreq, counts, counts.shape[0])
        info("Done!")
        # free some ram
        del counts

    info("Calculating entropy vector")
    popfreq_to_weights(popfreq, popfreq.shape[0], nsamples)
    weights = popfreq

    info("Writing weights to", outfile, end='... ')
    write_array(outfile, weights, name='weights')
    info("Done!")


def kernel_main():
    cli = '''
    USAGE:
        kwipy-kernel [options] OUTDIR WEIGHTFILE COUNTFILES ...

    OPTIONS:
        -c      Resume previous calculation to OUTDIR.
    '''
    opts = docopt(cli)
    outdir = opts['OUTDIR']
    weightfile = opts['WEIGHTFILE']
    countfiles = opts['COUNTFILES']
    resume = opts['-c']

    pairs = list(itl.combinations_with_replacement(countfiles, 2))

    # TODO: name kernel logs by the md5 of input filenames
    outfile = path.join(outdir, 'kernellog')
    if resume and path.exists(outfile):
        pairs_done = set()
        with open(outfile) as fh:
            for line in fh:
                a, b, kern = line.strip().split('\t')
                pairs_done.add((a, b))
        pairs = pairs.filter(lambda x: x not in pairs_done)
    else:
        with open(outfile, 'w') as fh:
            pass

    kernfn = partial(calc_kernel, weightfile)

    pool = Pool()
    for a, b, kernel in pool.imap_unordered(kernfn, pairs, 1):
        info(a, b, "=>",  kernel)
        with open(outfile, 'a') as kfh:
            print(a, b, kernel, sep='\t', file=kfh)


def distmat_main():
    cli = '''
    USAGE:
        kwipy-distmat [options] KERNELDIR

    OPTIONS:
        -k KERNEL   Output unnormalised kernel to KERNEL
        -d DIST     Output normalised distances to DIST
        -n NAME_RE  Regex to exract name from basename of count files. The
                    first match group is used. The regex must apply to all
                    files.  [default: (.*)\.(kct|bcz)]
        --keep      Keep kernel logs in KERNELDIR (which are normally removed)
    '''

    opts = docopt(cli)
    kerndir = opts['KERNELDIR']
    kernfile = opts['-k']
    distfile = opts['-d']
    namere = opts['-n']

    samples, kernmat = kernlog_to_kernmatrix(read_kernlog(kerndir), namere)
    check_psd(kernmat, "Kernel")
    distmat = kernel_to_distance(kernmat)
    check_psd(distmat, "Distance")

    if distfile:
        with open(distfile, 'w') as dfh:
            print_lsmat(distmat, samples, file=dfh)

    if kernfile:
        with open(kernfile, 'w') as dfh:
            print_lsmat(kernmat, samples, file=dfh)
    else:
        print_lsmat(kernmat, samples)
