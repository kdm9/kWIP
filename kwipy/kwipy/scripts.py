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
import bloscpack as bp
from docopt import docopt
import numpy as np
import numexpr as ne
from pymer import CountMinKmerCounter
import screed
from mpi4py import MPI

import itertools as itl
from sys import stderr, stdout
import os
from os import path

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

    popfreq = np.zeros((nt, ts), dtype=float)
    nsamples = len(countfiles)

    for countfile in countfiles:
        print("Loading",  countfile, end='... '); stdout.flush()
        counter = CountMinKmerCounter.read(countfile)
        counts = counter.array
        ne.evaluate('popfreq + where(counts > 0, 1, 0)', out=popfreq)
        print("Done!")

    print("Calculating entropy vector")

    ne.evaluate('popfreq / nsamples', out=popfreq)
    ne.evaluate('popfreq * log(popfreq) + (1 - popfreq) * log((1-popfreq))',
                out=popfreq)
    print("Writing to", outfile, end='... '); stdout.flush()
    bp.pack_ndarray_file(popfreq, outfile)
    print("Done!")

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

    popfreq = np.zeros((nt, ts), dtype=float)
    nsamples = len(countfiles)
    print(nsamples)

    for countfile in countfiles:
        print("Loading",  countfile, end='... '); stdout.flush()
        counter = CountMinKmerCounter.read(countfile)
        counts = counter.array
        ne.evaluate('popfreq + where(counts > 0, 1, 0)', out=popfreq)
        print("Done!")

    print("Calculating entropy vector")

    ne.evaluate('popfreq / nsamples', out=popfreq)
    ne.evaluate('-(popfreq * log(popfreq) + (1 - popfreq) * log((1-popfreq)))',
                out=popfreq)
    # workaround for numexpr's lack of isnan(): compare inequality, NaN != NaN
    ne.evaluate('where(popfreq != popfreq, 0, popfreq)', out=popfreq)
    print("Writing to", outfile, end='... '); stdout.flush()
    bp_args = bp.BloscArgs(cname='lz4', clevel=9, shuffle=False)
    bp.pack_ndarray_file(popfreq, outfile, blosc_args=bp_args)
    print("Done!")

def kernel_argparse():
    cli = '''

    USAGE:
        kwipy-kernelcalc OUTDIR WEIGHTFILE COUNTFILES ...
    '''

    opts = docopt(cli)
    outdir = opts['OUTDIR']
    weightfile = opts['WEIGHTFILE']
    countfiles = opts['COUNTFILES']
    return outdir, weightfile, countfiles


def kernel_mpi_main():
    outdir, weightfile, countfiles = kernel_argparse()
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    pairs = list(itl.combinations(countfiles, 2))

    if rank == 0:
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

    weights = bp.unpack_ndarray_file(weightfile)

    outfile = path.join(outdir, 'kernels_{}'.format(rank))
    with open(outfile, 'w') as fh:
        pass

    for afile, bfile in pairs:
        a = CountMinKmerCounter.read(afile).array
        b = CountMinKmerCounter.read(bfile).array

        print(rank, afile, bfile)
        kernel = ne.evaluate('sum(a * b * weights, axis=1)').min()

        with open(outfile, 'a') as kfh:
            print(afile, bfile, kernel, sep='\t', file=kfh)
