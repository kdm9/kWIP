from __future__ import print_function, division

'''Utilities for IO of blocked, compressed, HD5-encoded 1-D arrays'''
import numpy as np
import tables

from ._counter import KmerCounter
from .utils import chunked_fastq_iterator

'''The chunk size for Array IO'''
CHUNKSIZE = 2**20  # 1 MiB


class CounterSet(object):

    """CounterSet: A set of kmer counts for a number of samples"""

    def __init__(self, kmersize=21, cvsize=2e8, use_cms=False):
        """Create a CounterSet ready to accept samples

        :kmersize: int, length of kmers, aka $k$
        :cvsize: int, length of count vector
        :use_cms: bool, use a count-min sketch to record abundance

        `use_cms` enables the use of a CMS to minimise the impact of false
        positive lookups in the count vector. This requires significantly more
        memory.
        """
        self._kmersize = kmersize
        self._cvsize = cvsize
        self._use_cms = use_cms
        self._samples = {}
        self._current_sample = 0

    def _consume_file(self, filename):
        ctr = KmerCounter(self._kmersize, self._cvsize)
        for chunk in chunked_fastq_iterator(filename):
            ctr.consume(ctr)

    def count_sample(self, filename):
        pass


def write_sample(filename, array, name='array'):
    h5f = tables.open_file(filename, 'w')
    filt = tables.Filters(complib='blosc:blosclz', complevel=9, shuffle=True)
    h5f.create_earray('/', name, obj=array, shape=(0, FIXME),
                      chunkshape=(1, CHUNKSIZE,), filters=filt)
    h5f.close()


def read_array(filename, name='array'):
    h5f = tables.open_file(filename, 'r')
    arr = h5f.get_node('/' + name).read()
    h5f.close()
    return arr


def iter_blocks(filename, name='array'):
    h5f = tables.open_file(filename, 'r')
    arr = h5f.get_node('/' + name)
    for i in range(0, arr.shape[0], CHUNKSIZE):
        block = arr.read(i, i+CHUNKSIZE)
        yield i, block
    h5f.close()


def array_meta(filename, name='array'):
    h5f = tables.open_file(filename, 'r')
    arr = h5f.get_node('/' + name)
    meta = {
        'dtype': arr.dtype,
        'shape': arr.shape,
    }
    h5f.close()
    return meta

