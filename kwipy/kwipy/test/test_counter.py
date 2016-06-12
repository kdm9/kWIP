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

from __future__ import print_function, division

import pytest
from kwipy.counter import (
    iter_kmers,
    Counter
)

# de Bruijn DNA sequences of k={2,3}, i.e. contain all 2/3-mers once
K2_DBS = 'AACAGATCCGCTGGTTA'
K3_DBS = 'AAACAAGAATACCACGACTAGCAGGAGTATCATGATTCCCGCCTCGGCGTCTGCTTGGGTGTTTAA'


def do_test_iter_kmers(seq, expect, k):
    got = list(iter_kmers(seq, k))
    assert got == expect


def test_iter_kmers():
    '''Test iter_kmers(seq, k) with valid sequences'''
    # Valid
    seq = 'AACAGTA'
    expect = [0b0000, 0b0001, 0b0100, 0b0010, 0b1011, 0b1100]
    do_test_iter_kmers(seq, expect, 2)
    # Lower case
    seq = 'aacagta'
    do_test_iter_kmers(seq, expect, 2)
    # N splits
    seq = 'AACNGTA'
    expect = [0b0000, 0b0001,  # Skip 2 kmers with N
              0b1011, 0b1100]
    do_test_iter_kmers(seq, expect, 2)


def test_iter_kmers_bad():
    '''Test iter_kmers(seq, k) with bad sequences'''
    # Empty
    seq = ''
    do_test_iter_kmers(seq, [], 2)
    # < k
    seq = 'A'
    do_test_iter_kmers(seq, [], 2)
    # Too many Ns
    seq = 'ANANAN'
    do_test_iter_kmers(seq, [], 2)
    # All Ns
    seq = 'NNNNNN'
    do_test_iter_kmers(seq, [], 2)


def test_iter_kmers_err():
    '''Test iter_kmers(seq, k) with error cases'''
    # None
    seq = None
    with pytest.raises(TypeError):
        do_test_iter_kmers(seq, [], 2)
    # Bytes
    seq = b''
    with pytest.raises(TypeError):
        do_test_iter_kmers(seq, [], 2)


def test_counter_behaviour():
    '''Basic test of counting'''
    k = 3
    ctr = Counter(k=k, cvsize=1e5)
    ctr.consume(K3_DBS)
    for i in range(4**k):
        assert ctr.get(i) >= 1

    # The sum of the CV is not always the same as the number of k-mers. This is
    # because it is updated to the count-min sketch's estimate of the count of
    # an item.
    assert ctr.cv.sum() <= len(K3_DBS) - k + 1


def test_counter_nocms():
    '''Basic test of counting *without* a CMS'''
    k = 3
    ctr = Counter(k=k, cvsize=10, use_cms=False)
    ctr.consume(K3_DBS)

    # The sum of the CV **IS* the number of k-mers here, as the CMS is not used.
    assert ctr.cv.sum() == len(K3_DBS) - k + 1

    # Check that each k-mers is counted (at least) once
    for i in range(4**k):
        assert ctr.get(i) >= 1
