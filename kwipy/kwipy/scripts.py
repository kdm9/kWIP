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

import bcolz
from docopt import docopt
from mpi4py import MPI
import numpy as np
import numexpr as ne
from pymer import CountMinKmerCounter

import itertools as itl
from glob import glob
import os
from os import path
import re
from sys import stderr, stdout

from .kernelmath import *
from .logging import *
from .utils import (
    calc_weights,
    count_reads,
    print_lsmat,
)


def count_main():
    cli = '''
    USAGE:
        kwipy-count [options] OUTFILE READFILES ...

    OPTIONS:
        -k KSIZE    Kmer length [default: 20]
        -N NTAB     Number of tables [default: 4]
        -x TSIZE    Table size [default: 1e9]

    Counts k-mers in READFILES to a count-min sketch which is saved to OUTFILE.
    '''

    opts = docopt(cli)
    k = int(opts['-k'])
    nt = int(opts['-N'])
    ts = int(float(opts['-x']))
    outfile = opts['OUTFILE']
    readfiles = opts['READFILES']

    counts = count_reads(*readfiles, k=k, sketchshape=(nt, ts))

    info("Writing counts to", outfile)
    counts.write(outfile)


def weight_main():
    cli = '''
    USAGE:
        kwipy-weight [options] WEIGHTFILE COUNTFILES ...

    OPTIONS:
        -k KSIZE    Kmer length [default: 20]
        -N NTAB     Number of tables [default: 4]
        -x TSIZE    Table size [default: 1e9]
    '''

    opts = docopt(cli)
    k = int(opts['-k'])
    nt = int(opts['-N'])
    ts = int(float(opts['-x']))
    outfile = opts['WEIGHTFILE']
    countfiles = opts['COUNTFILES']

    weights = calc_weights(*countfiles, k=k, sketchshape=(nt, ts))

    info("Writing weights to", outfile, end='... ')
    weights = bcolz.carray(popfreq, rootdir=outfile, mode='w')
    weights.flush()
    info("Done!")


def kernel_mpi_main():
    cli = '''
    USAGE:
        kwipy-kernel-mpi [options] OUTDIR WEIGHTFILE COUNTFILES ...

    OPTIONS:
        -c      Resume previous calculation to OUTDIR

    '''
    opts = docopt(cli)
    outdir = opts['OUTDIR']
    weightfile = opts['WEIGHTFILE']
    countfiles = opts['COUNTFILES']
    resume = opts['-c']

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    pairs = list(itl.combinations_with_replacement(countfiles, 2))

    if rank == 0:
        if resume:
            existing_kernlogs = glob(path.join(outdir, 'kernellog_*'))
            pairs_done = set()
            for kl in existing_kernlogs:
                with open(kl) as fh:
                    for line in fh:
                        a, b, kern = line.strip().split('\t')
                        pairs_done.add((a, b))
            pairs = pairs.filter(lambda x: x not in pairs_done)
        size = comm.Get_size()
        if size > len(pairs):
            warn('Number of MPI ranks is greater than number of comparisons')
            warn('This is harmless but silly')
        pieces = [list() for x in range(size)]
        for i, pair in enumerate(pairs):
            pieces[i%size].append(pair)
    else:
        pieces = None

    pairs = comm.scatter(pieces, root=0)

    weights = bcolz.open(weightfile, mode='r')

    outfile = path.join(outdir, 'kernellog_{}'.format(rank))
    if not resume:
        with open(outfile, 'w') as fh:
            pass

    for afile, bfile in pairs:
        a = bcolz.open(afile, mode='r')
        b = bcolz.open(bfile, mode='r')

        print(rank, afile, bfile)
        kernel = bcolz.eval('sum(a * b * weights, axis=0)').min()

        with open(outfile, 'a') as kfh:
            print(afile, bfile, kernel, sep='\t', file=kfh)


def distmat_main():
    cli = '''
    USAGE:
        kwipy-distmat [options] KERNELDIR

    OPTIONS:
        -k KERNEL   Output unnormalised kernel to KERNEL
        -d DIST     Output normalised distances to DIST
        -n NAME_RE  Regex to exract name from basename of count files. The
                    first match group is used. The regex must apply to all
                    files.  [default: (.*)\.(blz|kpy)]
        --keep      Keep kernel logs in KERNELDIR (which are normally removed)
    '''

    opts = docopt(cli)
    outdir = opts['KERNELDIR']
    kernfile = opts['-k']
    distfile = opts['-d']

    def _file_to_name(filename):
        filename = path.basename(filename)
        try:
            return re.search(opts['-n'], filename).groups()[0]
        except (TypeError, AttributeError):
            warn("Name regex did not match '{}', using basename"
                 .format(filename))
            return path.splitext(filename)[0]

    kernlogs = glob(path.join(outdir, 'kernellog_*'))
    kernels = {}
    samples = set()
    for klfile in kernlogs:
        with open(klfile) as fh:
            for line in fh:
                a, b, kern = line.strip().split('\t')
                a, b = map(_file_to_name, (a, b))
                kern = float(kern)
                samples.update((a, b))  # add a & b
                try:
                    kernels[a][b] = kern
                except KeyError:
                    kernels[a] = {b: kern}

    # Make an ordered array
    num_samples = len(samples)
    kernmat = np.zeros((num_samples, num_samples), dtype=float)
    samples = list(sorted(samples))
    sample_idx = {s:i for i, s in enumerate(samples)}
    # square up the dictionary
    for a, b in itl.combinations_with_replacement(sorted(samples), 2):
        ai = sample_idx[a]
        bi = sample_idx[b]
        kernmat[ai, bi] = kernels[a][b]
        kernmat[bi, ai] = kernels[a][b]

    if is_psd(kernmat):
        info("Kernel matrix is positive semi-definite")
    else:
        warn("Kernel matrix is NOT positive semi-definite")

    if is_psd(distmat):
        info("Distance matrix is positive semi-definite")
    else:
        warn("Distance matrix is NOT positive semi-definite")

    distmat = kernel_to_distance(kernmat)
    if distfile:
        with open(distfile, 'w') as dfh:
            print_lsmat(distmat, samples, file=dfh)

    if kernfile:
        with open(kernfile, 'w') as dfh:
            print_lsmat(kernmat, samples, file=dfh)
    else:
        print_lsmat(kernmat, samples)
