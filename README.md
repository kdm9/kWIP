kmerclust
=========

Work on clustering samples using kmer analysis

Installation
============


You need `libkhmer.a` and the khmer C++ headers. This will be added to khmer as
part of [PR #788](https://github.com/ged-lab/khmer/pull/788), so will probably
make it into khmer 1.4.

Then:

    mkdir build && cd build
    cmake .. -DKHMER_ROOT=$HOME # or /usr/local
    make

