# kWIP

The k-mer Weighted Inner Product.

This software implements a *de novo*, *alignment free* measure of sample
genetic dissimilarity which operates upon raw sequencing reads. It is able to
calculate the genetic dissimilarity between samples without any reference
genome, and without assembling one.

[![Build Status](https://travis-ci.org/kdmurray91/kWIP.svg?branch=version2)](https://travis-ci.org/kdmurray91/kWIP)
[![Documentation Status](https://readthedocs.org/projects/kwip/badge/?version=latest)](https://kwip.readthedocs.org)

Full documentation is available at https://kwip.readthedocs.org

The code present in this branch is a python rewrite of the old kwip code. It
features better parallelism, better hashing code, and no dependence on khmer.
The docs and this readme might take a while to catch up with the changes I've
made -- please create a github issue to ask questions or if something seems
wrong.


# Installation


In a month or two when the code has stablised, I'll make conda packages for
kwipy and kmkm. For now, sorry, but you have to do the following. If you want
to save some time, you can do `conda env create -f kwip-conda-env.yml` to get
dependencies for kwipy/kmkm. Feel free to poke me, or make a pull request that
adds recipies youself if I'm being too slow!


```bash
pip3 install -e git+https://github.com/kmerkmer/kmkm.git#egg=kmkm
pip3 install -e git+https://github.com/kdmurray91/kwip.git@version2#egg=kwipy
```


# Usage

Something along the lines of:

```bash
kwipy-count -k 21 -v 1e9 -j $NCPUS samples.kc path/to/reads/*.fastq*
kwipy-kernel -D samples_distance.tsv -K samples_kernel.tsv -j $NCPUS samples.kc
```

Do `kwipy-count --help` or `kwipy-kernel --help` for more info. Again, actual
real documentation will happen at some point. Poke me if "some point" becomes
unhelfully far into the future. No question is too small or large.


# How it works

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


(this hasn't changed at all with version 2, althought results might differ very
slightly due to a different hashing system. In other words, your tree/PCA of
kWIP distances should be identical to the human eye, but kwip distances from
version 0.2.0 and version 2 might not be exactly the same.)


# License

 Public Licence V2. Note that some code from the old kWIP (somewhere in
./graveyard) is GPLv3+. If this is an issue, create a github issue.

# Publication

A publication describing kWIP has been [published in PLOS Computational
Biology](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1005727).
Please cite this, regardless of the kWIP version you use.
