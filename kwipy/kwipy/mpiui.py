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

from .counter import Counter
from .logging import (
    info,
    warn,
)
from .progress import ProgressLogger
from .utils import (
    calc_kernel,
    count_reads,
    parse_reads_with_precmd,
    stripext,
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
    with --precmd. The pre-command can be a shell pipeline combining the effects
    of multiple programs. Interleaved or single ended reads must be printed on
    stdout by the command(s). The pre-command uses the find/xargs/GNU Parallel
    convention of using a pair of '{}' to mark where the filename should be
    placed. Examples of a pre-command include:

        --precmd 'fastq-dump --split-spot --stdout {}'
        --precmd 'gzcat {} | trimit'
    '''

    opts = docopt(cli)
    k = int(opts['-k'])
    cvsize = int(float(opts['-v']))
    outdir = opts['OUTDIR']
    readfiles = opts['READFILES']
    use_cms = not opts['--no-cms']
    precmd = opts['--precmd']

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if rank == 0:
        size = comm.Get_size()
        if size > len(readfiles):
            warn('Number of MPI ranks is greater than number of comparisons')
            warn('This is harmless but silly')
        pieces = [list() for x in range(size)]
        for i, readfile in enumerate(readfiles):
            pieces[i % size].append(readfile)
    else:
        pieces = None

    our_readfiles = comm.scatter(pieces, root=0)

    for readfile in our_readfiles:
        base = stripext(readfile, ['fa', 'fq', 'fasta', 'fastq', 'gz', ])
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

    if rank == 0 and not path.isdir(outdir):
        mkdir(outdir)

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
            pieces[i % size].append(pair)
    else:
        pieces = None

    pairs = comm.scatter(pieces, root=0)

    outfile = path.join(outdir, 'kernellog_{}'.format(rank))
    if not resume:
        with open(outfile, 'w') as fh:
            pass

    for afile, bfile in pairs:
        a, b, k = calc_kernel(weightfile, (afile, bfile))
        print(a, b, k)
        with open(outfile, 'a') as kfh:
            print(a, b, k, sep='\t', file=kfh)
