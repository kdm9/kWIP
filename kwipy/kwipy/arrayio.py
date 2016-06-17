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


'''Utilities for IO of blocked, compressed, HD5-encoded 1-D arrays'''
import numpy as np
import tables


'''The chunk size for Array IO'''
CHUNKSIZE = 2**20  # 1 MiB


def write_array(filename, array, name='array'):
    h5f = tables.open_file(filename, 'w')
    filt = tables.Filters(complib='blosc:blosclz', complevel=9, shuffle=True)
    h5f.create_carray('/', name, obj=array, chunkshape=(CHUNKSIZE,),
                      filters=filt)
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
