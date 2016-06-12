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
