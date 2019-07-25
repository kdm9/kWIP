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
from fastqandfurious import fastqandfurious as fqf, _fastqandfurious as _fqf

from glob import glob
import itertools as itl
import os
from os import path
import re
import shlex
from shutil import rmtree
from subprocess import DEVNULL, PIPE, Popen

from .logging import (
    info,
    warn,
)
from .progress import ProgressLogger


READFILE_EXTS = ['gz', 'bz2', 'fa', 'fq', 'fasta', 'fastq', 'sra']


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
    mkdir(directory)


def chunked_fastq_iterator(filename, precmd=None):
    if precmd is not None:
        cmd = precmd.format(' "{}"'.format(filename))
        with Popen(cmd, shell=True, executable='/bin/bash', stdin=DEVNULL,
                   stdout=PIPE, stderr=None, universal_newlines=False) as proc:
            for read in FastxReader(proc.stdout):
                yield read
        if proc.returncode != 0:
            raise RuntimeError("Precommand exited with non-zero status.")
    else:
        with open(filename, 'rb') as fh:
            for read in FastxReader(fh):
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



