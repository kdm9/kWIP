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

from .logging import *


def print_lsmat(matrix, ids, file=None):
    print('', *ids, sep='\t', file=file)
    for rowidx, id in enumerate(ids):
        print(id, *list(matrix[rowidx,:]), sep='\t', file=file)


def count_reads(*readfiles, k=20, sketchshape=(4, 1e9)):
    counter = CountMinKmerCounter(k, sketchshape=sketchshape)

    for readfile in readfiles:
        info("Consuming:",  readfile)
        with screed.open(readfile) as reads:
            for i, read in enumerate(reads):
                counter.consume(str(read.sequence))
                if i % 10000 == 0:
                    progress(i/1000, 'K reads')
        progress(i/1000, 'K reads', end='\n')
    return counter


def calc_weights(*countfiles, k=20, sketchshape=(4, 1e9)):
    nt, ts = sketchshape
    popfreq = np.zeros((ts, nt), dtype=float)
    nsamples = len(countfiles)

    for countfile in countfiles:
        info("Loading",  countfile, end='... ')
        counts = bcolz.open(countfile, mode='r')
        # FIXME: the line below is due to a bug in bcolz, I think
        counts = np.array(counts).reshape(counts.shape)
        ne.evaluate('popfreq + where(counts > 0, 1, 0)', out=popfreq)
        info("Done!")

    info("Calculating entropy vector")

    # population counts to population freqs
    ne.evaluate('popfreq / nsamples', out=popfreq)
    # Population freqs to shannon entropy of pop freqs, i.e. weights
    ne.evaluate('-(popfreq * log(popfreq) + (1 - popfreq) * log((1-popfreq)))',
                out=popfreq)
    # workaround for numexpr's lack of isnan(): compare inequality, NaN != NaN
    ne.evaluate('where(popfreq != popfreq, 0, popfreq)', out=popfreq)
