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

from __future__ import print_function, division

import bcolz
import numpy as np
import numexpr as ne
from pymer import CountMinKmerCounter
import screed

from os import path

from .logging import (
    progress,
    info,
    warn,
)
from .counter import Counter
from .internals import (
    wipkernel,
    popfreq_add_sample,
    popfreq_to_weights
)


def print_lsmat(matrix, ids, file=None):
    print('', *ids, sep='\t', file=file)
    for rowidx, id in enumerate(ids):
        print(id, *list(matrix[rowidx, :]), sep='\t', file=file)


def count_reads(readfiles, k=20, cvsize=2e8):
    counter = Counter(k, cvsize)

    for readfile in readfiles:
        info("Consuming:",  readfile)
        with screed.open(readfile) as reads:
            for i, read in enumerate(reads):
                counter.consume(read.sequence)
                if i % 10000 == 0:
                    progress(i/1000, 'K reads')
        progress(i/1000, 'K reads', end='\n')
    return counter


def calc_weights(countfiles):
    popfreq = None
    nsamples = len(countfiles)

    for countfile in countfiles:
        info("Loading",  countfile, end='... ')
        counts = bcolz.open(countfile, mode='r')[:]
        if popfreq is None:
            popfreq = np.zeros(counts.shape, dtype=np.float32)
        popfreq_add_sample(popfreq, counts, counts.shape[0])
        info("Done!")

    info("Calculating entropy vector")
    popfreq_to_weights(popfreq, popfreq.shape[0], nsamples)
    return popfreq


def calc_kernel(weightsfile, ab):
    afile, bfile = ab
    abase = path.basename(afile)
    bbase = path.basename(bfile)

    A = bcolz.open(afile, mode='r')
    B = bcolz.open(bfile, mode='r')
    weights = bcolz.open(weightsfile, mode='r')

    kernel = 0.0

    bl = min([A.chunklen, B.chunklen, weights.chunklen])

    for a, b, w in zip(bcolz.iterblocks(A, blen=bl),
                       bcolz.iterblocks(B, blen=bl),
                       bcolz.iterblocks(weights, blen=bl)):
        kernel += wipkernel(a, b, w, bl)
    return abase, bbase, kernel
