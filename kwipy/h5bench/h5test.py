from __future__ import print_function, division
import h5py as h5
import numpy as np
import bcolz
import tables
from time import time
from sys import stdout

CHUNKSIZE = 1000000
NITER = 5
N = int(2e8)
DT = 'f4'

tables.set_blosc_max_threads(1)


def write_array_h5(filename, array, name='counts'):
    h5f = h5.File(filename, 'w')
    h5f.create_dataset(name, data=array, chunks=(CHUNKSIZE,), compression='lzf',
                       shuffle=False)
    h5f.close()


def write_array_blosc(filename, array, name='counts'):
    h5f = tables.open_file(filename, 'w')
    filt = tables.Filters(complib='blosc:lz4', complevel=1, shuffle=True)
    h5f.create_carray('/', name, obj=array, chunkshape=(CHUNKSIZE,),
                      filters=filt)
    h5f.close()


def write_array_lz4hc(filename, array, name='counts'):
    h5f = tables.open_file(filename, 'w')
    filt = tables.Filters(complib='blosc:lz4hc', complevel=9, shuffle=True)
    h5f.create_carray('/', name, obj=array, chunkshape=(CHUNKSIZE,),
                      filters=filt)
    h5f.close()


def write_array_h5s(filename, array, name='counts'):
    h5f = h5.File(filename, 'w')
    h5f.create_dataset(name, data=array, chunks=(CHUNKSIZE,), compression='lzf',
                       shuffle=True)
    h5f.close()


def write_array_bc(filename, array, name='counts'):
    arr = bcolz.carray(array=array, chunklen=CHUNKSIZE, rootdir=filename,
                       mode='w')
    arr.flush()


def read_array_h5(filename, _, name='counts'):
    h5f = h5.File(filename, 'w')
    arr = h5f[name][:]
    h5f.close()
    return arr


def read_array_pt(filename, _, name='counts'):
    h5f = tables.open_file(filename, 'r')
    arr = h5f.root.counts[:]
    h5f.close()
    return arr


def read_array_bc(filename, _, name='counts'):
    arr = bcolz.open(filename, mode='r')[:]
    return arr


def iter_blocks_bc(filename, _, name='counts'):
    f = bcolz.open(filename, mode='r')
    for block in bcolz.iterblocks(f):
        assert(np.sum(block) == CHUNKSIZE)


def iter_blocks_pt(filename, _, name='counts'):
    h5f = tables.open_file(filename, 'r')
    arr = h5f.root.counts
    csize = arr.chunkshape[0]
    for i in range(0, arr.shape[0], csize):
        assert(np.sum(arr[i:i+csize]) == csize)
    h5f.close()


def timeit(func, outf):
    array = np.ones(N, dtype=DT)
    start = time()
    for i in range(NITER):
        r = func(outf, array)
        print('.', end='')
        stdout.flush()
        if r is not None and i == 0:
            assert (array == r).all(), "Round trip failed"
    took = time() - start
    print()
    print(func.__name__, "(", outf, ") took", took, "sec, ", took/NITER,
          ' sec/iter')


# timeit(write_array_h5, "arr.h5")
# timeit(read_array_h5, "arr.h5")
# timeit(write_array_h5s, "arrs.h5")
# timeit(read_array_h5, "arrs.h5")
# timeit(write_array_lz4hc, "arrlz4hc.h5")
# timeit(read_array_pt, "arrlz4hc.h5")
timeit(write_array_blosc, "arrlz4.h5")
timeit(read_array_pt, "arrlz4.h5")
timeit(iter_blocks_pt, "arrlz4.h5")
timeit(write_array_bc, "arr.bcz")
timeit(read_array_bc, "arr.bcz")
timeit(iter_blocks_bc, "arr.bcz")
