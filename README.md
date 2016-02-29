kWIP
====

The k-mer Weighted Inner Product.

This software implements a *de novo*, *alignment free* measure of sample
genetic dissimilarity which operates upon raw sequencing reads. It is able to
calculate the genetic dissimilarity between samples without any reference
genome, and without assembling one.

[![Build Status](https://travis-ci.org/kdmurray91/kWIP.svg?branch=master)](https://travis-ci.org/kdmurray91/kWIP)
[![gplv3+](https://img.shields.io/badge/license-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl.html)
[![Documentation Status](https://readthedocs.org/projects/kwip/badge/?version=latest)](https://kwip.readthedocs.org)

Full documentation is available at https://kwip.readthedocs.org


Installation
============

### The full way

Dependencies:

- `zlib`
- `cmake>=2.8`
- Optionally, `liboxli` and `Eigen3` are required. These libraries are bundled
  with kWIP, and the internal copy will be used if system copies are not.
- A C++11 compiler that supports OpenMP (i.e. gcc >=4.8)

On Debian (or Debian derivatives) the dependencies of `kWIP` can be installed
with:

    sudo apt-get install zlib1g-dev cmake build-essential git

Then, to compile `kWIP`:

    git clone https://github.com/kdmurray91/kWIP.git
    cd kWIP
    mkdir build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=${HOME}
    make
    make test
    make install

The commands above assume you want to install kWIP to your home directory. This
is probably required on clusters, and necessary without root privileges. To
install to, e.g, `/usr/local/`, replace all occurrences of `$HOME` with your
preferred installation prefix.

### The easy way

Pre-compiled static binaries for 64-bit GNU/Linux are provided on the [GitHub
releases page](https://github.com/kdmurray91/kWIP/releases) in an archive named
`kwip-binaries_VERSION.tar.xz` or similar (where VERSION is the latest released
version). Please download this to your **GNU/Linux** machine, and:

    # Assuming you want to install to ~/bin/kwip
    PREFIX=$HOME
    cd $PREFIX

    # If ~/bin/ is not in $PATH, you won't be able to use kwip. Perform the
    # command below to ensure PATH is set correctly.
    echo "PATH=\"${PREFIX}/bin:\${PATH}\"" >> ~/.bashrc
    . ~/.bashrc

    # Below, repace all text in quotes with the URL to the archive on the
    # GitHub release page
    wget "URL FROM THE RELEASES PAGE"

    # Extract the precompile binaries and documentation
    tar xvf kwip-binaries*.tar.xz

    # Check your installation by typing
    kwip --help

    # The HTML documenation is available under ./share/doc/kwip
    ls $PREFIX/share/doc/kwip


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


License
=======

kWIP is Copyright 2015 Kevin Murray, and released under the GNU General Public
License version 3 (or any later version). See [THANKS.md](./THANKS.md) for a
list of contributors.

A publication describing kWIP is soon to be released publicly.
