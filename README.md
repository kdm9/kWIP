kmerclust
=========

Work on clustering samples using kmer analysis

Installation
============


You need `libkhmer.a` and the khmer C++ headers. These can be compiled from the
`master` branch of khmer:

    git clone https://github.com/ged-lab/khmer.git
    cd khmer/lib
    make install PREFIX=$HOME/

Then, to compile `kmerclust`:

    mkdir build && cd build
    cmake .. -DKHMER_ROOT=$HOME # or /usr/local
    make

