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

from blessings import Terminal
import bcolz
import bloscpack as bp
from docopt import docopt
import numpy as np
import numexpr as ne
from pymer import CountMinKmerCounter
import screed
from mpi4py import MPI

import itertools as itl
from glob import glob
import os
from os import path
import re
from sys import stderr, stdout


term = Terminal()

def progress(*args, file=stdout):
    file.write(term.move_x(0))
    file.write(term.clear_eol())
    print(*args, end='', file=file)
    file.flush()


def info(*args, file=stdout):
    print(term.blue, *args, term.normal, file=file)
    file.flush()


def warn(*args, file=stdout):
    print(term.bold_yellow, *args, term.normal, file=file)
    file.flush()


def hash_main():
    cli = '''

    USAGE:
        kwipy-hash [options] OUTFILE READFILES ...

    OPTIONS:
        -k KSIZE    Kmer length [default: 20]
        -N NTAB     Number of tables [default: 4]
        -x TSIZE    Table size [default: 1e9]
    '''

    opts = docopt(cli)
    k = int(opts['-k'])
    nt = int(opts['-N'])
    ts = int(float(opts['-x']))
    outfile = opts['OUTFILE']
    readfiles = opts['READFILES']

    counter = CountMinKmerCounter(k, sketchshape=(nt, ts))

    for readfile in readfiles:
        print("Consuming:",  readfile)
        with screed.open(readfile) as reads:
            for i, read in enumerate(reads):
                counter.consume(str(read.sequence))
                if i % 10000 == 0:
                    progress(i/1000, 'K reads')
        progress(i/1000, 'K reads')
        print()

    print("Writing to", outfile)
    counter.write(outfile)

def calcweight_main():
    cli = '''

    USAGE:
        kwipy-calcweight [options] WEIGHTFILE COUNTFILES ...

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

    popfreq = np.zeros((ts, nt), dtype=float)
    nsamples = len(countfiles)
    print(nsamples)

    for countfile in countfiles:
        print("Loading",  countfile, end='... '); stdout.flush()
        counts = bcolz.open(countfile, mode='r')
        counts = np.array(counts).reshape(counts.shape)
        ne.evaluate('popfreq + where(counts > 0, 1, 0)', out=popfreq)
        print("Done!")

    print("Calculating entropy vector")

    ne.evaluate('popfreq / nsamples', out=popfreq)
    ne.evaluate('-(popfreq * log(popfreq) + (1 - popfreq) * log((1-popfreq)))',
                out=popfreq)
    # workaround for numexpr's lack of isnan(): compare inequality, NaN != NaN
    ne.evaluate('where(popfreq != popfreq, 0, popfreq)', out=popfreq)
    print("Writing to", outfile, end='... '); stdout.flush()
    popfreq = bcolz.carray(popfreq, rootdir=outfile, mode='w')
    popfreq.flush()
    print("Done!")

def kernel_argparse():
    cli = '''

    USAGE:
        kwipy-kernelcalc [options] OUTDIR WEIGHTFILE COUNTFILES ...

    OPTIONS:
        -c      Resume previous calculation to OUTDIR
    '''

    opts = docopt(cli)
    outdir = opts['OUTDIR']
    weightfile = opts['WEIGHTFILE']
    countfiles = opts['COUNTFILES']
    resume = opts['-c']
    return outdir, weightfile, countfiles, resume


def kernel_mpi_main():
    outdir, weightfile, countfiles, resume = kernel_argparse()
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
        print(a.shape, b.shape, weights.shape)
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

    print(kernels)

    # square up the dictionary
    for a, b in itl.combinations_with_replacement(samples, 2):
        try:
            kernels[b][a] = kern
        except KeyError:
            kernels[b] = {a: kern}

    print(kernels)

    # Make an ordered array
    num_samples = len(samples)
    kernmat = np.zeros((num_samples, num_samples), dtype=float)
    for (ia, a), (ib, b) in zip(enumerate(sorted(samples)),
                                enumerate(sorted(samples))):
        kernmat[ia, ib] = kernels[a, b]

    print(kernmat)
