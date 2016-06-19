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


prof = True
if prof:
    import pstats
    import cProfile

    cProfile.runctx("do_it('/home/kevin/ws/seqs/paired/test_R1.fastq')",
                    globals(), locals(), "profile.prof")
    st = pstats.Stats("profile.prof")
    st.strip_dirs().sort_stats("cumtime").print_stats(1000)
else:
    if len(argv) != 2:
        print("USAGE: benchmark_counting.py <READFILE>", file=stderr)
        exit(1)
    do_it(argv[1])
if False:
    seq = 'ACGT'  # * 1000000
    counter = Counter(k=21, cvsize=1e9, use_cms=False)
    import line_profiler
    profile = line_profiler.LineProfiler(counter.consume)
    profile.runcall(counter.consume, seq)
    profile.print_stats()
