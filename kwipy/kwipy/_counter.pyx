import numpy as np
cimport numpy as np
cimport cython
from pymer._hash cimport iter_kmers


ctypedef unsigned long long int u64
ctypedef unsigned short u16
ctypedef unsigned char eltype
cpdef str dtype='u1'

cdef inline u64 mm64(u64 key, u64 seed):
    cdef u64 m = 0xc6a4a7935bd1e995
    cdef u64 r = 47

    cdef u64 h = seed ^ (8 * m)

    key *= m
    key ^= key >> r
    key *= m

    h ^= key
    h *= m

    return h


cdef class KmerCounter(object):
    cdef readonly u64 k
    cdef readonly u64 nt, ts, cvsize
    cdef u64 dtmax
    cdef readonly np.ndarray cv
    cdef np.ndarray cms
    cdef eltype[:] _cv
    cdef eltype[:, :] _cms

    def __init__(self, k, cvsize=2e8, use_cms=False):
        self.k = k
        if use_cms:
            self.nt, self.ts = (4, cvsize / 2)
        else:
            self.nt = self.ts = 0

        self.cvsize = cvsize
        self.dtmax = 2**16 - 1

        self.cv = np.zeros(int(cvsize), dtype=dtype)
        self.cms = np.zeros((self.nt, self.ts), dtype=dtype)
        self._cv = self.cv
        self._cms = self.cms

    @cython.nonecheck(False)
    @cython.boundscheck(False)
    @cython.overflowcheck(False)
    @cython.wraparound(False)
    cdef count(KmerCounter self, u64 item):
        cdef u64 count = 0xffffffffffffffff  # max value of u64
        cdef u64 hsh, cv_bin, v, i, tab

        cv_bin = mm64(item, 1) % self.cvsize
        if self.nt == 0:
            count = self._cv[cv_bin]
            if count < self.dtmax:
                count += 1
            self._cv[cv_bin] = count
            return count

        for tab in range(self.nt):
            hsh = mm64(item, tab + 1)
            i = hsh % self.ts
            v = self._cms[tab, i]
            if v < self.dtmax:
                v += 1
            if count > v:
                count = v
            self._cms[tab, i] = v
            self._cv[cv_bin] = count
        return count

    @cython.boundscheck(False)
    @cython.overflowcheck(False)
    @cython.wraparound(False)
    cpdef get(KmerCounter self, u64 item):
        cdef u64 hsh
        hsh = mm64(item, 1)
        return self._cv[hsh % self.cvsize]

    def consume(self, bytes seq not None):
        for kmer in iter_kmers(seq, self.k):
            self.count(kmer)

    def consume(self, list sequences not None):
        for seq in sequences:
            self.consume(seq)
