# Copyright (c) 2015-2019 Kevin Murray <foss@kdmurray.id.au>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import print_function, division

import numpy as np

from glob import glob
import itertools as itl
import os
from os import path
import re
import shlex
from shutil import rmtree
from sys import stdin, stdout, stderr
import functools
from multiprocessing import Pool

from .logging import (
    info,
    warn,
)


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



def print_lsmat(matrix, ids, file=stdin):
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


class Mapper(object):
    def __init__(self, njobs, chunksize=1):
        if njobs > 1:
            self.pool = Pool(njobs)
            self.mapper = functools.partial(self.pool.imap_unordered, chunksize=chunksize)
        else:
            self.mapper = map

    def __call__(self, *args, **kwargs):
        yield from self.mapper(*args, **kwargs)
