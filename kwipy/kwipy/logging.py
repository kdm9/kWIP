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

from sys import stdout, stderr

__all__ = ['progress', 'info', 'warn', 'error']


class MockTerm(object):
    def __getattr__(self, key):
        return ""


try:
    from blessings import Terminal
    term = Terminal()
except Exception as exc:
    term = MockTerm()


def progress(*args, file=stderr, end='\n'):
    file.write(term.cyan)
    print("  ... ", *args, end=end, file=file)
    file.write(term.normal)
    file.flush()


def info(*args, file=stderr, end='\n'):
    file.write(term.blue)
    print(*args, file=file, end=end)
    file.write(term.normal)
    file.flush()


def warn(*args, file=stderr, end='\n'):
    file.write(term.bold_yellow)
    print(*args, file=file, end=end)
    file.write(term.normal)
    file.flush()


def error(*args, file=stderr, end='\n'):
    file.write(term.bold_red)
    print(*args, file=file, end=end)
    file.write(term.normal)
    file.flush()
