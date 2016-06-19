#!/usr/bin/env python3
from __future__ import print_function, division
from kwipy.counter import Counter
from sys import argv, stdout, stderr, exit
from time import time
import screed


def do_it(fname):
    counter = Counter(k=21, cvsize=1e9, use_cms=False)
    with screed.open(fname) as reads:
        start = time()
        for i, read in enumerate(reads):
            counter.consume(read.sequence)
        took = time() - start
        print("Did", i + 1, "reads in {:0.2f}s".format(took))


if len(argv) == 2:
    do_it(argv[1])
else:
    import pstats
    import cProfile

    cProfile.runctx("do_it('/home/kevin/ws/seqs/paired/test_R1.fastq')",
                    globals(), locals(), "profile.prof")
    st = pstats.Stats("profile.prof")
    st.strip_dirs().sort_stats("cumtime").print_stats(1000)
