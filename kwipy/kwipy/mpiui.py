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
import bcolz

import itertools as itl
from glob import glob
from os import path, mkdir
import re
from sys import stderr, stdout
import shutil

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
    info,
    warn,
)
from .progress import ProgressLogger
from .utils import (
    calc_kernel,
    count_reads,
    parse_reads_with_precmd,
    read_kernlog,
    stripext,
    mpisplit,
    rmmkdir,
)


def count_mpi_main():
    cli = '''
    USAGE:
        kwipy-count-mpi [options] OUTDIR READFILES ...

    OPTIONS:
        -k KSIZE    Kmer length [default: 20]
        -v CVLEN    Count vector length [default: 1e9]
        -p PRECMD   Shell pipeline to run on input files before hashing.
        --no-cms    Disable the CMS counter, use only a count vector.
                    [default: False]

    Counts k-mers into individual count vectors, parallelised using MPI.

    Will use about 6 * CVLEN bytes of RAM per file (or 2x with --no-cms).

    An optional pre-counting command for e.g. QC or SRA dumping can be given
    with `-p`. The pre-command can be a shell pipeline combining the effects of
    multiple programs. Interleaved or single ended reads must be printed on
    stdout by the command(s). The pre-command uses the find/xargs/GNU Parallel
    convention of using a pair of '{}' to mark where the filename should be
    placed. Examples of a pre-command include:

        -p 'fastq-dump --split-spot --stdout {}'
        -p 'gzcat {} | trimit'
    '''

    opts = docopt(cli)
    k = int(opts['-k'])
    cvsize = int(float(opts['-v']))
    outdir = opts['OUTDIR']
    readfiles = opts['READFILES']
    use_cms = not opts['--no-cms']
    precmd = opts['--precmd']

    comm = MPI.COMM_WORLD
    readfiles = mpisplit(readfiles, comm)

    for readfile in readfiles:
        base = stripext(readfile, READFILE_EXTS)
        outfile = path.join(outdir, base + '.kct')
        counts = Counter(k=k, cvsize=cvsize, use_cms=use_cms)
        if precmd is None:
            count_reads(counts, [readfile, ])
        else:
            reads = parse_reads_with_precmd(readfile, precmd)
            for read in ProgressLogger(reads, 'reads', interval=10000):
                counts.consume(read.sequence)
        info("Writing counts to", outfile)
        counts.save(outfile)


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

    # In parallel, we summarise the counts into  a presence/absence table for
    # each rank
    popfreq = None
    for countfile in countfiles:
        info("Loading", countfile)
        counts = bcolz.open(countfile, mode='r')[:]
        if popfreq is None:
            popfreq = np.zeros(counts.shape, dtype=np.float32)
        popfreq_add_sample(popfreq, counts, counts.shape[0])

    temp_array = '{}_{}'.format(weightfile, rank)
    bcolz.carray(popfreq, rootdir=temp_array, mode='w',
                 chunklen=BCOLZ_CHUNKLEN).flush()

    # Then summarise the summarised counts and turn the freqs into weights
    if rank == 0:
        info("Summarising population counts")
        for i in range(1, comm.Get_size()):
            arrayfile = '{}_{}'.format(weightfile, i)
            theirpf = bcolz.open(arrayfile, mode='r')[:]
            inplace_append(popfreq, theirpf, popfreq.shape[0])
            shutil.rmtree(arrayfile)
        info("Calculating weights")
        popfreq_to_weights(popfreq, popfreq.shape[0], n_samples)
        bcolz.carray(popfreq, rootdir=weightfile, mode='w',
                     chunklen=BCOLZ_CHUNKLEN).flush()
        info("Wrote weights to", weightfile)
        shutil.rmtree('{}_{}'.format(weightfile, 0))


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

    pairs = list(itl.combinations_with_replacement(countfiles, 2))
    if resume:
        kernels = read_kernlog(outdir)
        pairs = pairs.filter(lambda p: p in kernels)
    else:
        rmmkdir(outdir)
    pairs = mpisplit(pairs, comm)

    outfile = path.join(outdir, 'kernellog_{}'.format(rank))

    for afile, bfile in pairs:
        a, b, k = calc_kernel(weightfile, (afile, bfile))
        if not quiet:
            print(a, b, k)
        with open(outfile, 'a') as kfh:
            print(a, b, k, sep='\t', file=kfh)
