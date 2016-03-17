import numpy as np
from numpy import linalg


def is_psd(matrix):
    '''Test if ``matrix`` is posititve semi-definite.'''
    v, w = linalg.eig(matrix)
    return (v > 0.0).all()


def normalise_kernel(K):
    '''Normalise a kernel matrix such that the diagonal is equal to 1.'''
    N = np.zeros_like(K)
    for i, j in itl.product(*map(range, K.shape)):
        N[i, j] = K[i, j] / np.sqrt(K[i,i] * K[j, j])
    return N


def kernel_to_distance(K):
    '''Convert a kernel matrix to a distance matrix.

    The kernel matrix need not be normalised, this is done internally.'''
    K = normalise_kern(kernel)
    D = np.zeros_like(kernel)
    for i, j in itl.product(*map(range, K.shape)):
        D[i,j] = (K[i,i] + K[j,j]) - 2*K[i,j]
    return D
