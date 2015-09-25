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

### The easy way

This will install `kWIP` to `~/bin`. If this is not what you want, use the full
install instructions below.

    git clone https://github.com/kdmurray91/kWIP.git
    cd kWIP
    ./util/simple-build.sh

This will build `khmer` and `kWIP`. Then, to install `kWIP`:

    cd build
    make install


### The full way

Dependencies:

- `liboxli.a` -- from the DIB-lab's [khmer](https://github.com/dib-lab/khmer)
- `cmake>=2.8`
- A C++11 compiler (gcc >=4.8, clang >=3.4).

`kWIP` depends upon `liboxli.a` and the khmer C++ headers during
compilation, version 2.0 or later.

    git clone https://github.com/dib-lab/khmer.git
    git checkout v2.0
    # Build the library and install it
    cd khmer/lib
    make install PREFIX=$HOME

Then, to compile `kWIP`:

    git clone https://github.com/kdmurray91/kWIP.git
    cd kWIP
    mkdir build && cd build  # Out-of-source build for sanity
    cmake .. -DOXLI_ROOT=$HOME
    make
    make test
    make install

The commands above assume you want to install kWIP and liboxli to your
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
