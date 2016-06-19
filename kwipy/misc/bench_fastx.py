#!/usr/bin/env python3
import screed
from kwipy.fastx import FastxReader, decompress
from sys import argv, exit
from time import time
import xxhash
from readfq import readfq
import io


if len(argv) != 2:
    print("USAGE: bench_fastx <READFILE>")
    exit(1)
filename = argv[1]


def count_seqlen(reads):
    totlen = 0
    xxh = xxhash.xxh64()
    for read in reads:
        totlen += len(read.sequence)
        xxh.update(read.sequence)
    return totlen, xxh.hexdigest()

start = time()
with open(filename, 'rb') as fh:
    reads = FastxReader(fh)
    cnt = count_seqlen(reads)
    print("kwipy:", cnt, "took {:0.2f}s".format(time() - start))

start = time()
with open(filename, 'rb') as fh:
    totlen = 0
    xxh = xxhash.xxh64()
    for n, s, q in readfq(io.TextIOWrapper(decompress(fh))):
        totlen += len(s)
        xxh.update(s)
    cnt = (totlen, xxh.hexdigest())
    print("readfq:", cnt, "took {:0.2f}s".format(time() - start))

start = time()
with screed.open(filename) as reads:
    cnt = count_seqlen(reads)
    print("screed:", cnt, "took {:0.2f}s".format(time() - start))


print('performing sanity check')
with screed.open(filename) as s_reads, open(filename, 'rb') as fh:
    k_reads = FastxReader(fh)
    for s, k in zip(s_reads, k_reads):
        assert bytes(s.name, 'ascii') == k.name
        assert bytes(s.sequence, 'ascii') == k.sequence
        assert bytes(s.quality, 'ascii') == k.quality
print('passed!')
