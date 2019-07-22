# Copyright 2016 Kevin Murray <kdmfoss@gmail.com>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.


import numpy as np
cimport numpy as cnp
cimport cython


ctypedef unsigned long long int u64
ctypedef unsigned int u32


@cython.boundscheck(False)
def hash_to_kmer(int h, int k):
    '''Convert a hash value at given k to the string representation.
    '''
    cdef const char *nts = 'ACGT'
    cdef u64 nt
    kmer = []
    for x in range(k):
        nt = (h >> (2*x)) & 0x03
        kmer.append(chr(nts[nt]))
    return ''.join(reversed(kmer))


cdef rc_kmer(u64 x, int k):
    # Reverse
    x = (x & 0x3333333333333333) <<  2 | (x & 0xCCCCCCCCCCCCCCCC) >>  2
    x = (x & 0x0F0F0F0F0F0F0F0F) <<  4 | (x & 0xF0F0F0F0F0F0F0F0) >>  4
    x = (x & 0x00FF00FF00FF00FF) <<  8 | (x & 0xFF00FF00FF00FF00) >>  8
    x = (x & 0x0000FFFF0000FFFF) << 16 | (x & 0xFFFF0000FFFF0000) >> 16
    x = (x & 0x00000000FFFFFFFF) << 32 | (x & 0xFFFFFFFF00000000) >> 32
    # Complement
    x = ~x
    # Shift back within bitmask
    x >>= 64 - 2 * k
    return x

def iter_kmers(str seq, int k, bint canonical=False, bint mask=False):
    '''Iterator over hashed k-mers in a string DNA sequence.
    '''
    cdef u64 n
    cdef u64 bitmask = 2**(2*k)-1  # Set lowest 2*k bits
    cdef u64 h = 0

    # For each kmer's end nucleotide, bit-shift, add the end and yield
    cdef u64 skip = 0
    for end in range(len(seq)):
        nt = seq[end]
        if skip > 0:
            skip -= 1
        if nt == 'A' or (nt == 'a' and not mask):
            n = 0
        elif nt == 'C' or (nt == 'c' and not mask):
            n = 1
        elif nt == 'G' or (nt == 'g' and not mask):
            n = 2
        elif nt == 'T' or (nt == 't' and not mask):
            n = 3
        else:
            skip = k
            continue
        h = ((h << 2) | n) & bitmask
        if end >= k - 1 and skip == 0:
            # Only yield once an entire kmer has been loaded into h
            if canonical:
                yield min(h, rc_kmer(h, k))
            else:
                yield h
