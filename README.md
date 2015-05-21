kmerclust
=========

Work on clustering samples using kmer analysis

[![Build Status](http://biojenkins.anu.edu.au/job/kmerclust/badge/icon)](http://biojenkins.anu.edu.au/job/kmerclust/)
[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

Installation
============


You need `libkhmer.a` and the Khmer C++ headers. These can be compiled from
version 1.4 or greater of Khmer:

    git clone https://github.com/dib-lab/khmer.git
    git checkout v1.4
    cd khmer/lib
    make install PREFIX=$HOME/

Then, to compile `kmerclust`:

    mkdir build && cd build
    cmake .. -DKHMER_ROOT=$HOME
    make

The commands above assume you want to install kmerclust and libkhmer to your
home directory. This is probably required on clusters, and necessary without
root privileges. To install to, e.g, `/usr/local/`, replace all occurrences of
`$HOME` with your preferred installation prefix.
