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

See the [``kWIP``
Documentation](http://kwip.readthedocs.org/en/latest/installation.html)

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
License version 3 (or any later version).

Publication
===========

A publication describing kWIP has been [published in PLOS Computational Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005727)
