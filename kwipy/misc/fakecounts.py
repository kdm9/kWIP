#!/usr/bin/env python3
import bcolz
import numpy as np

import sys


arrays = [
    [0, 0, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 1],
    [0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 1, 1, 1],
    [0, 0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1],
    [0, 1, 0, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1],
    [0, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1],
    [0, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 0, 1, 1],
    [0, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0],
]

outdir = sys.argv[1]

for i, arr in enumerate(arrays):
    arr = np.array(arr, dtype=np.uint16)

    outf = "{}/{}.kct".format(outdir, i)
    bcolz.carray(arr, rootdir=outf, mode='w').flush()
    print(i, arr)
