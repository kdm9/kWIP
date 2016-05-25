import numpy as np
cimport numpy as np
from bcolz import carray
from pymer._hash import iter_kmers
cimport cython
from .constants import BCOLZ_CHUNKLEN


ctypedef unsigned long long int u64


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


cdef class Counter(object):
    cdef readonly u64 k
    cdef readonly u64 nt, ts, cvsize
    cdef u64 dtmax
    cdef np.ndarray cms, cv

    def __init__(self, k, cvsize=2e8):
        self.k = k
        self.nt, self.ts = (4, cvsize / 2)

        self.cvsize = cvsize
        dtype='u2'
        self.dtmax = 2**16 - 1

        if self.nt < 1 or self.nt > 10:
            raise ValueError("Too many or few tables. must be 1 <= nt <= 10")

        self.cms = np.zeros((self.nt, self.ts), dtype=dtype)
        self.cv = np.zeros(int(cvsize), dtype=dtype)

    @cython.boundscheck(False)
    @cython.overflowcheck(False)
    @cython.wraparound(False)
    cdef count(Counter self, u64 item):
        cdef u64 minv = <u64>-1
        cdef u64 hsh, v, i

        for tab in range(self.nt):
            hsh = mm64(item, tab)
            i = hsh % self.ts
            v = self.cms[tab, i]
            if v <= self.dtmax:
                v += 1
            if minv > v:
                minv = v
            self.cms[tab, i] = v
        self.cv[hsh % self.cvsize] = minv
        return minv

    def consume(self, str seq not None):
        cdef long long hashv = 0
        for kmer in iter_kmers(seq, self.k):
            self.count(kmer)

    def save(self, str filename not None):
        carray(self.cv, rootdir=filename, mode='w',
               chunklen=BCOLZ_CHUNKLEN).flush()
