#!/usr/bin/env python3
from __future__ import print_function, division
from collections import defaultdict
import pandas as pd
import numpy as np

from sys import argv, stderr, stdout


def parse_log(logfile):
    dists = defaultdict(dict)

    with open(logfile) as fh:
        for line in fh:
            line = line.rstrip('\r\n')
            if line.startswith("#") or not line:
                continue
            s1, s2, idx, dist = line.split()
            dists[s1][s2] = float(dist)
            dists[s2][s1] = float(dist)
    return dists


def print_matrix(dists, file=stdout):
    samples = list(sorted(dists.keys()))
    names = [''] + samples
    print(*names, sep='\t', file=file)

    for sample1 in samples:
        line = [sample1]
        for sample2 in samples:
            line.append(dists[sample1][sample2])
        print(*line, sep='\t', file=file)


def main():
    logfile = argv[1]
    dists = parse_log(logfile)
    print_matrix(dists)

if __name__ == '__main__':
    main()
