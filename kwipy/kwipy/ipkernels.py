from __future__ import print_function, division, absolute_import

import numpy as np
from numpy import linalg
import itertools as itl
import zarr
import functools


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
    f = ((X>1).sum(axis=0).astype(float) /  n)
    #import numexpr as ne
    #h = ne.evaluate("-(f * log2(f) + (1-f) * log2(1-f))")
    h = np.nan_to_num(-(f * np.log2(f) + (1-f) * np.log2(1-f)))
    #K = X @ np.diag(h) @ X.T
    W = X * h
    K = W @ W.T
    return K


def ipkernel(X):
    return X @ X.T


def chunked_calc(idx, arrayfile, kernel):
    start, end = idx
    arr = zarr.open(arrayfile)
    X = arr[:, start:end].astype(float)
    K =  kernel(X)
    return K


def eval_kernel(array, kernel, maxp=None, mapper=map, chunksize=None):
    X = zarr.open(array)
    c = X.chunks[1] if chunksize is None else chunksize
    n, p = X.shape
    if maxp is None or maxp > p:
        maxp = p
    indicies = [(s, min(maxp, s+c)) for s in range(0, maxp, c)]
    f = functools.partial(chunked_calc, arrayfile=array, kernel=kernel)
    K = None
    for k in mapper(f, indicies):
        if K is None:
            K = k
            continue
        K += k
    return K
