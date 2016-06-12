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

from __future__ import print_function, division, absolute_import

import numpy as np
from numpy import linalg
import itertools as itl


def is_psd(matrix):
    '''Test if ``matrix`` is posititve semi-definite.'''
    v, w = linalg.eig(matrix)
    return (v >= 0.0).all()


def normalise_kernel(K):
    '''Normalise a kernel matrix such that the diagonal is equal to 1.'''
    N = np.zeros_like(K, dtype=float)
    for i, j in itl.product(*map(range, K.shape)):
        N[i, j] = K[i, j] / np.sqrt(K[i, i] * K[j, j])
    return N


def kernel_to_distance(K):
    '''Convert a kernel matrix to a distance matrix.

    The kernel matrix need not be normalised, this is done internally.'''
    if not is_psd(K):
        raise UserWarning("Converting non-euclidean kernel matrix to distance")
    K = normalise_kernel(K)
    D = np.zeros_like(K, dtype=float)
    for i, j in itl.product(*map(range, K.shape)):
        D[i, j] = (K[i, i] + K[j, j]) - 2*K[i, j]
    return D
