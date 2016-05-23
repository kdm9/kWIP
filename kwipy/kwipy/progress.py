# Copyright 2016 Kevin Murray <spam@kdmurray.id.au>
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:
#
# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
# BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
# ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

from __future__ import print_function, division, absolute_import
import progressbar
from progressbar import UnknownLength, ProgressBar, FormatLabel
import sys
from sys import stdout, stderr
from time import time, sleep
from random import random
from collections import Iterable, Iterator

__all__ = ['ProgressLogger', ]


def human_readable(count, maxsuffix='M'):
    '''Returns a human-readable representation of the size of something'''
    for unit in ['', 'K', 'M', 'G', 'T', 'P', 'E', 'Z', 'Y']:
        if abs(count) < 1000.0 or unit == maxsuffix:
            break
        count /= 1000.0
    if unit != '':
        return "{i:.1f}{unit}".format(i=count, unit=unit)
    return "{i:.0f}{unit}".format(i=count, unit=unit)


class ProgressLogger(object):
    '''A logger of progress.

    Example:
    --------

        log = ProgressLogger(interval=100)
        for x in range(100):
            sleep(0.1)
            log.update(x)
        log.done()

    Useful as a wrapper around an iterable:

        for x in ProgressLogger(range(100), interval=10):
            sleep(0.1)
    '''
    moving_chars = "-\\|/"

    def __init__(self, iterable=None, interval=1000, human_readable=True,
                 poll_interval=None, file=stderr, format='{value} ({speed}/s)',
                 noun=''):
        self.interval = interval
        self.human_readable = human_readable
        self.poll_interval = poll_interval
        self.start_time = None
        self.last_update = None
        self.file = file
        self.fmt = format
        self.last = None  # holds last call to update
        self.noun = noun

        if self.file.isatty():
            self.print_prefix = '\r\x1b[K\x1b[36m --> '  # clear line, start cyan
            self.print_suffix = '\x1b[m\x0f'  # normal text, \r
        else:
            self.print_prefix = '   ... '
            self.print_suffix = '\n'

        if iterable is not None:
            if not isinstance(iterable, Iterable):
                raise ValueError("Non-iterable given as iterable kwarg")
            if not isinstance(iterable, Iterator):
                iterable = iter(iterable)
        self.iterable = iterable
        self.i = None

    def __iter__(self):
        if self.iterable is not None:
            self.i = 0
            return self
        else:
            return None

    def __next__(self):
        if self.iterable is None:
            raise ValueError("No iterable given")
        try:
            self.update(self.i, self.noun)
            v = next(self.iterable)
            self.i += 1
            return v
        except StopIteration as exc:
            self.done()
            raise exc

    def _should_update(self, i):
        if self.start_time is None:
            # first call to update()
            return True
        if self.interval is not None and self.interval > 0:
            if i % self.interval == 0:
                return True
        if self.poll_interval is not None and self.last_update is not None:
            if time() - self.last_update > self.poll_interval:
                return True
        return False

    def _update(self, i, args):
        now = time()
        elapsed = now - self.start_time

        value = i
        speed = i / elapsed
        if self.human_readable:
            value = human_readable(i)
            speed = human_readable(i/elapsed)

        updatestr = self.fmt.format(value=value, time=elapsed, speed=speed)

        self.file.write(self.print_prefix)
        print(updatestr, *args, end='', file=self.file)
        self.file.write(self.print_suffix)
        self.file.flush()
        self.last_update = now

    def update(self, i, *args):
        if not self.start_time:
            # First update
            if i == 0:
                self.first = 0
            else:
                self.first = 1
            self.start_time = time()
        self.last = (i, args)
        if not self._should_update(i):
            return
        self._update(i, args)

    def done(self):
        i, args = self.last
        self._update(i, args)
        print(file=self.file)
