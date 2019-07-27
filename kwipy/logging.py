# Copyright (c) 2015-2019 Kevin Murray <foss@kdmurray.id.au>
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
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
    print("    -", *args, end=end, file=file)
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
