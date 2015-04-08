kmerclust
=========

Work on clustering samples using kmer analysis

[![Build Status](http://biojenkins.anu.edu.au/job/kmerclust/badge/icon)](http://biojenkins.anu.edu.au/job/kmerclust/)
[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

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

