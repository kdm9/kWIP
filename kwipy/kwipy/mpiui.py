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

from mpi4py import MPI
from docopt import docopt
import numpy as np

import itertools as itl
from glob import glob
from os import path, mkdir
import re
from sys import stderr, stdout
import shutil

from . import cliargs
from .counter import Counter
from .constants import (
    READFILE_EXTS,
    BCOLZ_CHUNKLEN,
)
from .internals import (
    popfreq_add_sample,
    popfreq_to_weights,
    inplace_append,
)
from .logging import (
    error,
    info,
    warn,
)
from .progress import ProgressLogger
from .utils import (
    calc_kernel,
    count_reads,
    needs_update,
    read_array,
    read_kernlog,
    rmmkdir,
    stripext,
    write_array,
)


def mpisplit(things, comm):
    '''Split `things` into a chunk for each MPI rank.'''
    pieces = None
    rank = comm.Get_rank()
    if rank == 0:
        size = comm.Get_size()
        if size > len(things):
            warn('Number of MPI ranks is greater than number of items')
            warn('This is harmless but silly')
        pieces = [list() for x in range(size)]
        for i, thing in enumerate(things):
            pieces[i % size].append(thing)
    return comm.scatter(pieces, root=0)


def count_mpi_main():
    parser = cliargs.count_args()
    args = parser.parse_args()

    comm = MPI.COMM_WORLD
    readfiles = mpisplit(args.readfiles, comm)

    for readfile in readfiles:
        base = stripext(readfile, READFILE_EXTS)
        outfile = path.join(args.outfile, base + '.kct')
        if not needs_update(readfile, outfile):
            info("Skipping", readfile, "as it is older than", outfile)
            continue
        info("Consuming reads from", readfile)
        counts = Counter(k=args.ksize, cvsize=int(args.cvsize),
                         use_cms=args.use_cms)
        try:
            count_reads(counts, [readfile, ], precmd=args.precmd)
            info("Writing counts to", outfile)
            counts.save(outfile)
        except Exception:
            break


def weight_mpi_main():
    cli = '''
    USAGE:
        kwipy-weight-mpi WEIGHTFILE COUNTFILES ...
    '''
    opts = docopt(cli)
    weightfile = opts['WEIGHTFILE']
    countfiles = opts['COUNTFILES']
    n_samples = len(countfiles)

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    countfiles = mpisplit(countfiles, comm)

    # In parallel, we summarise the counts into a presence/absence table for
    # each rank
    popfreq = None
    for countfile in countfiles:
        info("Loading", countfile)
        counts = read_array(countfile)
        if popfreq is None:
            popfreq = np.zeros(counts.shape, dtype=np.float32)
        popfreq_add_sample(popfreq, counts, counts.shape[0])

    temp_array = '{}_{}'.format(weightfile, rank)
    write_array(temp_array, popfreq)

    comm.Barrier()

    # Then summarise the summarised counts and turn the freqs into weights
    if rank == 0:
        info("Summarising population counts")
        for i in range(1, comm.Get_size()):
            arrayfile = '{}_{}'.format(weightfile, i)
            theirpf = read_array(arrayfile)
            inplace_append(popfreq, theirpf, popfreq.shape[0])
        info("Calculating weights")
        popfreq_to_weights(popfreq, popfreq.shape[0], n_samples)
        write_array(weightfile, popfreq)
        info("Wrote weights to", weightfile)
        for i in range(comm.Get_size()):
            shutil.rmtree('{}_{}'.format(weightfile, i))


def kernel_mpi_main():
    cli = '''
    USAGE:
        kwipy-kernel-mpi [options] OUTDIR WEIGHTFILE COUNTFILES ...

    OPTIONS:
        -Q      Quiet mode [default: false]
        -c      Resume previous calculation to OUTDIR
    '''
    opts = docopt(cli)
    outdir = opts['OUTDIR']
    weightfile = opts['WEIGHTFILE']
    quiet = opts['-Q']
    countfiles = opts['COUNTFILES']
    resume = opts['-c']

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0 and not resume:
        rmmkdir(outdir)
    if rank == 0 and not path.isdir(outdir):
        mkdir(outdir)

    comm.Barrier()

    pairs = list(itl.combinations_with_replacement(countfiles, 2))
    if resume:
        kernels = read_kernlog(outdir)
        pairs = pairs.filter(lambda p: p in kernels)
    pairs = mpisplit(pairs, comm)

    outfile = path.join(outdir, 'kernellog_{}'.format(rank))

    comm.Barrier()

    for afile, bfile in pairs:
        a, b, k = calc_kernel(weightfile, (afile, bfile))
        if not quiet:
            print(a, b, k)
        with open(outfile, 'a') as kfh:
            print(a, b, k, sep='\t', file=kfh)
