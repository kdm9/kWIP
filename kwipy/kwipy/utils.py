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

from glob import glob
import itertools as itl
from os import path, mkdir
import re
from shutil import rmtree

from .logging import (
    info,
    warn,
)
from .progress import ProgressLogger
from .counter import Counter
from .kernelmath import is_psd
from .internals import (
    wipkernel,
    popfreq_add_sample,
    popfreq_to_weights
)


def stripext(filename, extensions):
    filename = path.basename(filename)
    for ext in extensions:
        # Should handle cases like stripext(.tar.gz, [.tar, .tar.gz]).
        if filename == ext:
            break
        # Add a leading dot to the extension if it doesn't have one
        if not ext.startswith("."):
            ext = "." + ext
        # strip away extension
        if filename.endswith(ext):
            filename = filename[:-len(ext)]
    return filename


def check_psd(matrix, name):
    if is_psd(matrix):
        info(name, "matrix is positive semi-definite")
    else:
        warn(name, "matrix is NOT positive semi-definite")


def file_to_name(filename, namere=r'(.*)(\.[^.]+)?'):
    filename = path.basename(filename)
    try:
        return re.search(namere, filename).groups()[0]
    except (TypeError, AttributeError):
        warn("Name regex did not match", filename, "using basename")
        return filename


def print_lsmat(matrix, ids, file=None):
    print('', *ids, sep='\t', file=file)
    for rowidx, id in enumerate(ids):
        print(id, *list(matrix[rowidx, :]), sep='\t', file=file)


def rmmkdir(directory):
    if path.exists(directory):
        rmtree(directory)
    mkdir(directory)


def count_reads(readfiles, k=20, cvsize=2e8):
    counter = Counter(k, cvsize)

    for readfile in readfiles:
        info("Consuming:",  readfile)
        with screed.open(readfile) as reads:
            for read in ProgressLogger(reads, 'reads', interval=10000):
                counter.consume(read.sequence)
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
        # free some ram
        del counts

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


def read_kernlog(outdir):
    kernlines = []
    kernlogs = glob(path.join(outdir, 'kernellog*'))
    for kfn in kernlogs:
        with open(kfn) as fh:
            kernlines.extend(fh.readlines())
    return kernlines


def kernlog_to_kernmatrix(kernlog_lines, namere):
    kernels = {}
    samples = set()
    for line in kernlog_lines:
        a, b, kern = line.rstrip().split('\t')
        a = file_to_name(a, namere)
        b = file_to_name(b, namere)
        samples.add(a)
        samples.add(b)
        kern = float(kern)
        try:
            kernels[a][b] = kern
        except KeyError:
            kernels[a] = {b: kern}

    # Make an ordered array
    num_samples = len(samples)
    kernmat = np.zeros((num_samples, num_samples), dtype=float)
    samples = list(sorted(samples))
    sample_idx = {s: i for i, s in enumerate(samples)}
    # square up the dictionary
    for a, b in itl.combinations_with_replacement(samples, 2):
        ai = sample_idx[a]
        bi = sample_idx[b]
        kernmat[ai, bi] = kernels[a][b]
        kernmat[bi, ai] = kernels[a][b]
    return samples, kernmat
