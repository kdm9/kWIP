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

import numpy as np
import screed

from glob import glob
import itertools as itl
import os
from os import path
import re
import shlex
from shutil import rmtree
from subprocess import DEVNULL, PIPE, Popen

from .arrayio import (
    read_array,
    write_array,
    iter_blocks,
)
from .internals import (
    wipkernel,
)
from .kernelmath import is_psd
from .logging import (
    info,
    warn,
)
from .progress import ProgressLogger


def stripext(filename, extensions):
    """Strips each extension in ``extensions`` from ``filename``

    Args
    ----
    filename: str
        path (possibly with leading directories) to a file
    extensions: list of str
        extensions to strip. Extensions may start with a '.'; if they do not,
        one will be added.

    Returns
    -------
    str
        The basename of ``filename`` with any extensions in ``extensions``
        removed.
    """
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


def needs_update(source, dest):
    '''Checks if file ``source`` is newer than ``dest``.

    Args
    ----
    source: str
        path to a source file.
    dest: list of str
        File, derived from ``source``, that may be older than ``source``

    Returns
    -------
    bool:
        returns mtime(source) > mtime(dest) if both source and dest exist,
        otherwise True
    '''
    if path.exists(source) and path.exists(dest):
        return path.getmtime(source) > path.getmtime(dest)
    else:
        return True


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


def mkdir(directory):
    try:
        os.mkdir(directory)
    except OSError:
        pass


def rmmkdir(directory):
    try:
        if path.exists(directory):
            rmtree(directory)
    except OSError:
        pass
    os.mkdir(directory)


def parse_reads(filename, precmd=None):
    if precmd is not None:
        cmd = precmd.format(' "{}"'.format(filename))
        with Popen(cmd, shell=True, executable='/bin/bash', stdin=DEVNULL,
                   stdout=PIPE, stderr=None, universal_newlines=True) as proc:
            for read in screed.fastq.fastq_iter(proc.stdout):
                yield read
        if proc.returncode != 0:
            raise RuntimeError("Precommand exited with non-zero status.")
    else:
        with screed.open(filename) as reads:
            for read in reads:
                yield read


def count_reads(counter, readfiles, precmd=None, interval=50000):
    for readfile in readfiles:
        try:
            reads = parse_reads(readfile, precmd=precmd)
            for read in ProgressLogger(reads, 'reads', interval=interval):
                counter.consume(read.sequence)
        except IOError as exc:
            raise RuntimeError("An error occured while reading from", readfile,
                               ". The error was:\n\t", str(exc))


def calc_kernel(weightsfile, ab):
    afile, bfile = ab
    abase = path.basename(afile)
    bbase = path.basename(bfile)

    A = iter_blocks(afile, name='counts')
    B = iter_blocks(bfile, name='counts')
    W = iter_blocks(weightsfile, name='weights')

    kernel = 0.0

    # A/B/W all are generators returning a tuple of (i, block) where i is the
    # index in the big array of block[0].
    # Hence, ai == a's i, aa = a's array and so on.
    for (ai, aa), (bi, ba), (wi, wa) in zip(A, B, W):
        assert ai == bi == wi
        kernel += wipkernel(aa, ba, wa, aa.shape[0])
    return abase, bbase, kernel


def read_kernlog(outdir):
    kernlines = []
    kernlogs = glob(path.join(outdir, 'kernellog*'))
    for kfn in kernlogs:
        with open(kfn) as fh:
            kernlines.extend(fh.readlines())
    kernels = {}
    for line in kernlines:
        a, b, k = line.strip().split()
        kernels[(a, b)] = float(k)
    return kernels


def kernlog_to_kernmatrix(kernels_raw, namere):
    kernels = {}
    samples = set()
    for (a, b), kern in kernels_raw.items():
        a = file_to_name(a, namere)
        b = file_to_name(b, namere)
        samples.add(a)
        samples.add(b)
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
