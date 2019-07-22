
import numpy as np
from numpy import linalg
import itertools as itl



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
