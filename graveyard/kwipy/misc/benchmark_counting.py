#!/usr/bin/env python3
from __future__ import print_function, division
from kwipy.counter import Counter
from sys import argv, stdout, stderr, exit
from time import time
import screed

profile = True
seq = 'ACGT' * 25


def do_it():
    counter = Counter(k=21, cvsize=1e9, use_cms=False)
    start = time()
    for i in range(500000):
        counter.consume(seq)
    done_c = time()
    counter.save("test.h5")
    done_s = time()
    print("Count: {:0.2f}s, Save: {:0.2f}s".format(done_c - start,
                                                   done_s - done_c))


if not profile:
    do_it()
else:
    import pstats
    import cProfile

    cProfile.runctx("do_it()", globals(), locals(), "profile.prof")
    st = pstats.Stats("profile.prof")
    st.strip_dirs().sort_stats("cumtime").print_stats(30)
