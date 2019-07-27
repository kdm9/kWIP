# Copyright (c) 2015-2019 Kevin Murray <foss@kdmurray.id.au>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
from __future__ import print_function, division, absolute_import

import numpy as np
from numpy import linalg
import zarr

import functools
import itertools as itl
from multiprocessing import Pool
import warnings


def is_psd(matrix):
    '''Test if ``matrix`` is posititve semi-definite.'''
    v, w = linalg.eig(matrix)
    return (v >= -1e-6).all()


def normalise_kernel(K):
    '''Normalise a kernel matrix such that the diagonal is equal to 1.'''
    N = np.zeros_like(K, dtype=float)
    for i, j in itl.product(*map(range, K.shape)):
        N[i, j] = K[i, j] / np.sqrt(K[i, i] * K[j, j])
    return N


def kernel_to_distance(K):
    '''Convert a kernel matrix to a distance matrix.

    The kernel matrix need not be normalised, this is done internally.'''
    if not is_psd(K):
        raise UserWarning("Converting non-euclidean kernel matrix to distance")
    K = normalise_kernel(K)
    D = np.zeros_like(K, dtype=float)
    for i, j in itl.product(*map(range, K.shape)):
        D[i, j] = (K[i, i] + K[j, j]) - 2*K[i, j]
    return D


def wipkernel(X):
    n, p = X.shape
    f = ((X>0).sum(axis=0).astype(float) /  n)
    np.seterr(all="ignore")
    h = np.nan_to_num(-(f * np.log2(f) + (1-f) * np.log2(1-f)))                                                                                                             
    W = X * h
    K = W @ W.T
    return K

def ipkernel(X):
    return X @ X.T


def chunked_calc(idx, array, kernel):
    start, end = idx
    X = array[:, start:end]
    K =  kernel(X)
    return K


def eval_kernel(array, kernel, maxp=None, jobs=1, chunksize=None):
    c = array.chunks[1] if chunksize is None else chunksize
    n, p = array.shape
    if maxp is None or maxp > p:
        maxp = p
    indicies = [(s, min(maxp, s+c)) for s in range(0, maxp, c)]
    f = functools.partial(chunked_calc, array=array, kernel=kernel)
    K = None
    if jobs > 1:
        pool = Pool(jobs)
        mapper = functools.partial(pool.imap_unordered, chunksize=1)
    else:
        mapper = map
    for k in mapper(f, indicies):
        if K is None:
            K = k
            continue
        K += k
    return K
