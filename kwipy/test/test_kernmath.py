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
import numpy as np
from numpy.testing import assert_allclose, assert_equal

from kwipy.kernelmath import (
    is_psd,
    normalise_kernel,
    kernel_to_distance,
)

KERN = np.array([[8, 4, 0],
                 [4, 8, 0],
                 [0, 0, 8]])

DIST = np.array([[0, 1, 2],
                 [1, 0, 2],
                 [2, 2, 0]], dtype=float)

NOTPSD = np.array([[6, 1, 2],
                   [8, 0, 3],
                   [6, 2, 4]])


def test_is_psd():
    assert is_psd(KERN)
    assert not is_psd(NOTPSD)


def test_kernel_to_distance():
    '''Incorporates tests of normalise_kernel'''
    assert_allclose(kernel_to_distance(KERN), DIST)
    with pytest.raises(UserWarning):
        kernel_to_distance(NOTPSD)
