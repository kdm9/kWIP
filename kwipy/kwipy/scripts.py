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

from docopt import docopt
import numpy as np

import itertools as itl
from glob import glob
import os
from os import path
from sys import stderr, stdout
from multiprocessing import Pool
from functools import partial

from .arrayio import (
    read_array,
    write_array,
    iter_blocks,
)
from . import cliargs
from .counter import Counter
from .kernelmath import (
    normalise_kernel,
    kernel_to_distance,
)
from .internals import (
    popfreq_add_sample,
    popfreq_to_weights
)
from .logging import (
    info,
    warn,
)
from .utils import (
    calc_kernel,
    check_psd,
    count_reads,
    file_to_name,
    kernlog_to_kernmatrix,
    mkdir,
    read_kernlog,
    print_lsmat,
)


def count_main():
    parser = cliargs.count_args()
    args = parser.parse_args()

    mkdir(path.dirname(args.outfile))
    counter = Counter(k=args.ksize, cvsize=args.cvsize, use_cms=args.use_cms)
    try:
        count_reads(counter, args.readfiles, precmd=args.precmd)
        info("Writing counts to", args.outfile)
        counter.save(args.outfile)
    except Exception:
        pass


def weight_main():
    cli = '''
    USAGE:
        kwipy-weight WEIGHTFILE COUNTFILES ...
    '''

    opts = docopt(cli)
    outfile = opts['WEIGHTFILE']
    countfiles = opts['COUNTFILES']

    popfreq = None
    nsamples = len(countfiles)

    for countfile in countfiles:
        info("Loading",  countfile, end='... ')
        counts = read_array(countfiles, name='counts')
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
