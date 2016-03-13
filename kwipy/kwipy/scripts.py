from blessings import Terminal
from docopt import docopt
from pymer import CountMinKmerCounter
import screed
from sys import stderr, stdout


term = Terminal()

def progress(*args, file=stdout):
    print(term.move_x(0), term.clear_eol(), *args, end='', file=file)
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
        print()

    print("Writing to", outfile)
    counter.write(outfile)


