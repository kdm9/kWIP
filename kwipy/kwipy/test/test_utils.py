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

import numpy as np
from numpy.testing import assert_allclose, assert_equal
from tempdir import run_in_tempdir

from kwipy import utils


def test_stripext():
    assert utils.stripext('test.tar.gz', ['tar', 'tar.gz']) == "test"
    assert utils.stripext('path/to/test.tar.gz', ['.tar', '.tar.gz']) == "test"
    assert utils.stripext('.tar.gz', ['.tar', '.tar.gz']) == ".tar.gz"

def test_print_lsmat(capfd):
    utils.print_lsmat(np.array([[1, 0], [0, 1]]), ids=["A", "B"])
    out, err = capfd.readouterr()
    assert err == ""
    assert out == "\tA\tB\nA\t1\t0\nB\t0\t1\n"


@run_in_tempdir()
def test_array_io():
    arr = np.arange(1000, dtype=float)
    utils.write_array("test.bcz", arr)
    new_arr = utils.read_array("test.bcz")
    assert_equal(arr, new_arr)
