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

from blessings import Terminal
import bloscpack as bp
from docopt import docopt
import numpy as np
import numexpr as ne
from pymer import CountMinKmerCounter
import screed
from sys import stderr, stdout

term = Terminal()

def progress(*args, file=stdout):
    file.write(term.move_x(0))
    file.write(term.clear_eol())
    print(*args, end='', file=file)
    file.flush()


def info(*args, file=stdout):
    print(term.blue, *args, term.normal, file=file)
    file.flush()


def warn(*args, file=stdout):
    print(term.bold_yellow, *args, term.normal, file=file)
    file.flush()


def hash_main():
    cli = '''

    USAGE:
        kwipy-hash [options] OUTFILE READFILES ...

    OPTIONS:
        -k KSIZE    Kmer length [default: 20]
        -N NTAB     Number of tables [default: 4]
        -x TSIZE    Table size [default: 1e9]
    '''

    opts = docopt(cli)
    k = int(opts['-k'])
    nt = int(opts['-N'])
    ts = int(float(opts['-x']))
    outfile = opts['OUTFILE']
    readfiles = opts['READFILES']

    counter = CountMinKmerCounter(k, sketchshape=(nt, ts))

    for readfile in readfiles:
        print("Consuming:",  readfile)
        with screed.open(readfile) as reads:
            for i, read in enumerate(reads):
                counter.consume(str(read.sequence))
                if i % 10000 == 0:
                    progress(i/1000, 'K reads')
        progress(i/1000, 'K reads')
        print()

    print("Writing to", outfile)
    counter.write(outfile)


def calcweight_main():
    cli = '''

    USAGE:
        kwipy-calcweight [options] WEIGHTFILE COUNTFILES ...

    OPTIONS:
        -k KSIZE    Kmer length [default: 20]
        -N NTAB     Number of tables [default: 4]
        -x TSIZE    Table size [default: 1e9]
    '''

    opts = docopt(cli)
    k = int(opts['-k'])
    nt = int(opts['-N'])
    ts = int(float(opts['-x']))
    outfile = opts['WEIGHTFILE']
    countfiles = opts['COUNTFILES']

    popfreq = np.zeros((nt, ts), dtype=float)
    nsamples = len(countfiles)

    for countfile in countfiles:
        print("Loading",  countfile, end='... '); stdout.flush()
        counter = CountMinKmerCounter.read(countfile)
        counts = counter.array
        ne.evaluate('popfreq + where(counts > 0, 1, 0)', out=popfreq)
        print("Done!")

    print("Calculating entropy vector")

    ne.evaluate('popfreq / nsamples', out=popfreq)
    ne.evaluate('popfreq * log(popfreq) + (1 - popfreq) * log((1-popfreq))',
                out=popfreq)
    print("Writing to", outfile, end='... '); stdout.flush()
    bp.pack_ndarray_file(popfreq, outfile)
    print("Done!")
