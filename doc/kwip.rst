========
``kWIP``
========

The :math:`k`-mer Weighted Inner Product


Overview
~~~~~~~~

``kWIP`` is a method for calculating genetic similarity between samples. Unlike
similar alternatives, e.g. SNP-based distance calculation, ``kWIP`` operates
directly upon next-gen sequencing reads. ``kWIP`` works by decomposing
sequencing reads to short :math:`k`-mers, hashing these :math:`k`-mers and
performing pairwise distance calculation between these sample :math:`k`-mer
hashes. We use khmer from the DIB lab, UC Davis to hash sequencing reads.
``kWIP`` calculates the distance between samples in a computationally efficient
manner, and generates a distance matrix which may be used by downstream tools.
The power of ``kWIP`` comes from the weighting applied across different hash
values, which decreases the effect of erroneous, rare or over-abundant
:math:`k`-mers while focusing on :math:`k`-mers which give the most insight
into the similarity of samples.


Installation
~~~~~~~~~~~~

Dependencies:

- liboxli
- cmake (>=2.8)
- A C++11 compiler (e.g., gcc >=4.8).
- A fully POSIX-compliant operating system. ``kWIP`` is fully tested only on
  Debian GNU/Linux, however should work on any modern GNU/Linux and even Apple
  Mac OS X.

kWIP depends upon liboxli from the `khmer project
<https://github.com/dib-lab/khmer>`_, version 2.0 or greater.  In Debian
GNU/Linux, it should soon be possible to install ``liboxli-dev`` using ``sudo
apt-get install liboxli-dev``. On all other systems, the following commands
will compile and install ``liboxli``.

::

    git clone https://github.com/dib-lab/khmer.git
    cd khmer/lib
    make install PREFIX=$HOME

Then, to compile kWIP:

::

    git clone https://github.com/kdmurray91/kWIP.git
    cd kWIP/build
    cmake .. -DKHMER_ROOT=$HOME
    make
    make test
    make install

The commands above assume you want to install kWIP and libkhmer to your home
directory. This is probably required on clusters, and necessary without root
privileges. To install to, e.g, ``/usr/local/``, replace all occurrences of
``$HOME`` with your preferred installation prefix.


``kWIP`` CLI Usage
~~~~~~~~~~~~~~~~~~

