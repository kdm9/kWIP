import numpy as np
cimport numpy as np
cimport cython

ctypedef unsigned short u16
ctypedef unsigned long long u64

@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.wraparound(False)
def popfreq_add_sample(np.ndarray[float, ndim=1] popfreq, np.ndarray[u16, ndim=1]  sample, u64 nelem):
    for i in range(nelem):
        if sample[i] > 0:
            popfreq[i] += 1

@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.wraparound(False)
def popfreq_to_weights(np.ndarray[float, ndim=1] popfreq, u64 nelem, float nsamples):
    """In-place convertion of population abundance counts to entropy weights."""
    for i in range(nelem):
        freq = popfreq[i] / nsamples
        if freq == 0. or freq == 1.:
            popfreq[i] = 0.
            continue
        freq1m = 1 - freq
        ent = -((freq * np.log2(freq)) + (freq1m) * np.log2(freq1m))
        popfreq[i] = ent


@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.wraparound(False)
def wipkernel(np.ndarray[u16, ndim=1] A, np.ndarray[u16, ndim=1] B,
              np.ndarray[float, ndim=1] weights, u64 nelem):
    """Weighted inner product kernel"""
    cdef double kernel = 0.0
    for i in range(nelem):
        kernel += A[i] * B[i] * weights[i]
    return kernel
