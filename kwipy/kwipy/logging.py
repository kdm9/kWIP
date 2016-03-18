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

from blessings import Terminal

from sys import stdout, stderr


term = Terminal()


def progress(*args, file=stdout, end=''):
    file.write(term.move_x(0))
    file.write(term.clear_eol())
    file.write(term.cyan)
    print(*args, end=end, file=file)
    file.write(term.normal)
    file.flush()


def info(*args, file=stdout, end='\n'):
    file.write(term.blue)
    print(*args, file=file, end=end)
    file.write(term.normal)
    file.flush()


def warn(*args, file=stderr, end='\n'):
    file.write(term.bold_yellow)
    print(*args, file=file, end=end)
    file.write(term.normal)
    file.flush()
