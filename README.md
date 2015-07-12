kWIP
====

The k-mer Weighted Inner Product.

This software implements a *de novo*, *alignment free* measure of sample
genetic dissimilarity which operates upon raw sequencing reads. It is able to
calculate the genetic dissimilarity between samples without any reference
genome, and without assembling one.

[![Build Status](http://biojenkins.anu.edu.au/job/kwip/badge/icon)](http://biojenkins.anu.edu.au/job/kwip/)
[![Build Status](https://travis-ci.org/kdmurray91/kWIP.svg?branch=master)](https://travis-ci.org/kdmurray91/kWIP)
[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)
[![Documentation Status](https://readthedocs.org/projects/kwip/badge/?version=latest)](https://kwip.readthedocs.org)

Full documentation is available at https://kwip.readthedocs.org


Installation
============

Dependencies:

- `libkhmer.a`
- `cmake>=2.8`
- A C++11 compiler (gcc >=4.8, clang >=3.4).

`kWIP` depends upon `libkhmer.a` and uses the khmer C++ headers during
compilation. These should be compiled from the `master` branch of khmer, as
bugs relevant to our use case exist in the latest released version (v1.4.1).

    git clone https://github.com/dib-lab/khmer.git
    # Build the library and install it
    cd khmer/lib
    make install PREFIX=$HOME

Then, to compile `kWIP`:

    git clone https://github.com/kdmurray91/kWIP.git
    cd kWIP
    mkdir build && cd build  # Out-of-source build for sanity
    cmake .. -DKHMER_ROOT=$HOME
    make
    make test
    make install

The commands above assume you want to install kWIP and libkhmer to your
home directory. This is probably required on clusters, and necessary without
root privileges. To install to, e.g, `/usr/local/`, replace all occurrences of
`$HOME` with your preferred installation prefix.


How it works
============

kWIP works by decomposing sequencing reads to short
[k-mers](https://en.wikipedia.org/wiki/K-mer),
[hashing](https://en.wikipedia.org/wiki/Hash_function) these k-mers and
performing pairwise distance calculation between these sample k-mer hashes. We
use [`khmer`](https://github.com/dib-lab/khmer) from the DIB lab, UC Davis to
hash sequencing reads. `KWIP` calculates the distance between samples in a
computationally efficient manner, and generates a distance matrix which may be
used by downstream tools. The power of `kWIP` comes from the weighting applied
across different hash values, which decreases the effect of erroneous, rare or
over-abundant k-mers while focusing on k-mers which give the most insight into
the similarity of samples.

More detail will be added here soon.


License
=======

kWIP is Copyright 2015 Kevin Murray, and released under the GNU General Public
License version 3 (or any later version). See [THANKS.md](./THANKS.md) for a
list of contributors.

A publication describing kWIP is soon to be released publicly.


Miscellany
==========

- **Naming**: The working title of this software was `kmerclust`, any
  references to `kmerclust` are a result of my poor `sed` skills.
