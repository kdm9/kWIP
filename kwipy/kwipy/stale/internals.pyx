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

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport log2, isnan

ctypedef unsigned short u16
ctypedef unsigned long long u64


@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.wraparound(False)
def popfreq_add_sample(np.ndarray[float, ndim=1] popfreq,
                       np.ndarray[u16, ndim=1] sample, u64 nelem):
    for i in range(nelem):
        if sample[i] > 0:
            popfreq[i] += 1


@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.wraparound(False)
def inplace_append(np.ndarray[float, ndim=1] A, np.ndarray[float, ndim=1] B,
                   u64 nelem):
    for i in range(nelem):
        A[i] += B[i]


@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def popfreq_to_weights(np.ndarray[float, ndim=1] popfreq, u64 nelem,
                       float nsamples):
    """In-place convertion of population abundance counts to entropy weights."""
    cdef float freq, freq1m, ent
    for i in range(nelem):
        freq = popfreq[i] / nsamples
        if freq == 0. or freq == 1.:
            # entropy is zero, but is caclulated as NAN
            popfreq[i] = 0.
            continue
        freq1m = 1 - freq
        ent = -((freq * log2(freq)) + (freq1m) * log2(freq1m))
        if isnan(ent):
            print(freq, freq1m, ent)
        popfreq[i] = ent


@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.wraparound(False)
def wipkernel(np.ndarray[u16, ndim=1] A, np.ndarray[u16, ndim=1] B,
              np.ndarray[float, ndim=1] weights, u64 nelem):
    """Weighted inner product kernel"""
    cdef double kernel = 0.0
    cdef double a, b, w
    for i in range(nelem):
        # These are separate here to ensure fpmath is used in kernel calculation
        a = A[i]
        b = B[i]
        w = weights[i]
        kernel += a * b * w
    return kernel


@cython.boundscheck(False)
@cython.overflowcheck(False)
@cython.wraparound(False)
def ipkernel(np.ndarray[u16, ndim=1] A, np.ndarray[u16, ndim=1] B, u64 nelem):
    """Inner product kernel"""
    cdef double kernel = 0.0
    cdef double a, b
    for i in range(nelem):
        # These are separate here to ensure fpmath is used in kernel calculation
        a = A[i]
        b = B[i]
        kernel += a * b
    return kernel
