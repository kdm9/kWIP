========
``kWIP``
========

The :math:`k`-mer Weighted Inner Product


Overview
--------

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
------------

Dependencies:

- cmake (>=2.8)
- A C++11 compiler that supports OpenMP (e.g., g++ >=4.8).
- A fully POSIX-compliant operating system. ``kWIP`` is fully tested only on
  Debian GNU/Linux, however should work on any modern GNU/Linux and even Apple
  Mac OS X.

On Debian (or Debian derivatives) the dependencies of `kWIP` can be installed
with:

::
    sudo apt-get install zlib1g cmake build-essential

Then, to compile `kWIP`:

::
    git clone https://github.com/kdmurray91/kWIP.git
    cd kWIP
    mkdir build && cd build
    cmake ..
    make
    make test
    make install

The commands above assume you want to install kWIP to your home directory. This
is probably required on clusters, and necessary without root privileges. To
install elsewhere (e.g ``/usr/local/``), replace all occurrences of ``$HOME``
with your preferred installation prefix.


``kWIP`` CLI Usage
------------------

::

    USAGE: kwip [options] sample1 sample2 ... sampleN

    OPTIONS:
    -t, --threads       Number of threads to utilise. [default N_CPUS]
    -k, --kernel        Output file for the kernel matrix. [default None]
    -d, --distance      Output file for the distance matrix. [default stdout]
    -U, --unweighted    Use the unweighted inner proudct kernel. [default off]
    -w, --weights       Bin weight vector file (input, or output w/ -C).
    -C, --calc-weights  Calculate only the bin weight vector, not kernel matrix.
    -h, --help          Print this help message.
    -V, --version       Print the version string.
    -v, --verbose       Increase verbosity. May or may not acutally do anything.
    -q, --quiet         Execute silently but for errors.


The ``kwip`` executable is the core of ``kWIP``; its help statement is
reproduced above. This program operates on the saved Countgraphs of ``khmer``.
One can run with or without the entropy weighting, using the ``-U`` parameter
to disable weighting.

An example command could be:

::

    kwip \
        -t 4 \                        # Use 4 threads
        -k rice.kern \                # Output kernel matrix to ./rice.kern
        -d rice.dist \                # Output distance matrix to ./rice.dist
        ./hashes/rice_sample_*.ct.gz  # Path to sample hashes, with wildcard

Note that this is purely illustrative and won't run as-is due to the in-line
comments. Were it to run, it would calculate the Weighted Innner Product (WIP)
kernel pairwise between all samples given as arguments, utilising four threads
and saving the raw kernel matrix to rice.kern and the normalised distance
matrix to rice.dist.


The Concepts Behind ``kWIP``
----------------------------

The inner product between two vectors is directly related to the distance
between the vectors in Euclidean space. This has been utilised several times in
Bioinformatics to implement measures of genetic similarity between two
sequences, including the :math:`D2` statistic. Traditionally, the software
which implement these and similar algorithms operate on known genetic
sequences, e.g. those taken from a reference genome. ``kWIP``'s innovation is
to weight the inner product operation by a weight vector, and to derive weights
in a way which minimises the noise inherent in next-gen sequencing datasets
while maximising the signal of genetic distance between samples.
