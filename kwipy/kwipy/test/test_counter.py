import pytest
from kwipy.counter import (
    iter_kmers,
)


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


def test_counter():
    pass
