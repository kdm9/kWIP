
import numpy as np
cimport numpy as np
cimport cython
from .arrayio import write_array


ctypedef unsigned long long int u64
ctypedef unsigned short u16


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


def iter_kmers(str seq not None, int k):
    '''Iterator over hashed k-mers in a string DNA sequence.
    '''
    cdef u64 n
    cdef u64 bitmask = (1 << (2 * k)) - 1  # 4**k - 1
    cdef u64 h = 0

    # For each end nucleotide, bit-shift left, OR w/ the end nuc and yield
    cdef u64 skip = 0
    for end in range(len(seq)):
        nt = seq[end]
        if skip > 0:
            skip -= 1
        if nt == 'A' or nt == 'a':
            n = 0
        elif nt == 'C' or nt == 'c':
            n = 1
        elif nt == 'G' or nt == 'g':
            n = 2
        elif nt == 'T' or nt == 't':
            n = 3
        else:
            skip = k
            continue
        h = ((h << 2) | n) & bitmask
        if end >= k - 1 and skip == 0:
            # Only yield once an entire kmer has been loaded into h
            yield h


cdef class Counter(object):
    cdef readonly u64 k
    cdef readonly u64 nt, ts, cvsize
    cdef u64 dtmax
    cdef readonly np.ndarray cv
    cdef np.ndarray cms
    cdef u16[:] _cv
    cdef u16[:,:] _cms

    def __init__(self, k, cvsize=2e8, use_cms=True):
        self.k = k
        if use_cms:
            self.nt, self.ts = (4, cvsize / 2)
        else:
            self.nt = self.ts = 0

        self.cvsize = cvsize
        dtype='u2'
        self.dtmax = 2**16 - 1

        self.cv = np.zeros(int(cvsize), dtype=dtype)
        self.cms = np.zeros((self.nt, self.ts), dtype=dtype)
        self._cv = self.cv
        self._cms = self.cms

    @cython.boundscheck(False)
    @cython.overflowcheck(False)
    @cython.wraparound(False)
    cdef count(Counter self, u64 item):
        cdef u64 count = 0xffffffffffffffff  # max value of u64
        cdef u64 hsh, cv_bin, v, i

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
    cpdef get(Counter self, u64 item):
        cdef u64 hsh
        hsh = mm64(item, 1)
        return self._cv[hsh % self.cvsize]

    def consume(self, str seq not None):
        for kmer in iter_kmers(seq, self.k):
            self.count(kmer)

    def save(self, str filename not None):
        write_array(filename, self.cv, name='counts')
