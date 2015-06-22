kmerclust
=========

Work on clustering samples using kmer analysis

[![Build Status](http://biojenkins.anu.edu.au/job/kmerclust/badge/icon)](http://biojenkins.anu.edu.au/job/kmerclust/)
[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)

Installation
============


You need `libkhmer.a` and the khmer C++ headers. These should be compiled from
the `master` branch of khmer, as bugs relevant to our use case exist in the
latest released version v1.4.1.

    git clone https://github.com/dib-lab/khmer.git
    # Build the library and install it
    cd khmer/lib
    make install PREFIX=$HOME

Then, to compile `kmerclust`:

    git clone https://github.com/kdmurray91/kmerclust.git
    cd kmerclust
    mkdir build && cd build  # Out-of-source build for sanity
    cmake .. -DKHMER_ROOT=$HOME
    make
    make test
    make install

The commands above assume you want to install kmerclust and libkhmer to your
home directory. This is probably required on clusters, and necessary without
root privileges. To install to, e.g, `/usr/local/`, replace all occurrences of
`$HOME` with your preferred installation prefix.
