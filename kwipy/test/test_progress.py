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

from kwipy.progress import (
    human_readable,
    ProgressLogger,
)

def test_human_readable():
    units = ['', 'K', 'M']

    for i, u in enumerate(units):
        i = 1000 ** i
        assert human_readable(i) == "1{}".format(u)

    assert human_readable(1e9, maxsuffix='M') == '1000M'
    assert human_readable(1e9, maxsuffix='G') == '1G'
