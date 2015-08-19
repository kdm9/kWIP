===================================
Example ``kWIP`` Analysis Protocols
===================================


*Oryza sativa* grouping
-----------------------

This example uses data from the `3000 Rice Genomes Project
<http://http//www.gigasciencejournal.com/content/3/1/7>`_. The rice genome
isn't tiny, so you'll want a machine with at least 8gb RAM, and a least 2
cores.


Requirements
^^^^^^^^^^^^

- virutalenvwrapper
- git
- the SRA toolkit
- ``kWIP``, installed such that the ``kwip`` binary is in ``$PATH``
- We also need the ``img.R`` script which comes with ``kWIP``, but isn't
  installed. Hereafter I assume that you have this script in the analysis
  directory, ``~/rice_kwip``. Copy it with ``cp util/img.R ~/rice_kwip`` from
  the directory you ``git clone``-d or unziped ``kWIP`` into.


Setup
^^^^^

Create a working directory

::

    mkdir ~/rice_kwip
    cd ~/rice_kwip

Install ``SRApy`` and ``khmer`` using pip. We will require the master branch of
khmer until the 2.0 release, as bugs exist in the latest release uploaded to
PyPI.

::

    mkvirtualenv rice-kwip
    pip install srapy
    pip install -e git+https://github.com/dib-lab/khmer.git

Create the file ``sra_run_ids.txt`` with the following contents

::

    ERR626208
    ERR626209
    ERR626210
    ERR626211
    ERR619511
    ERR619512
    ERR619513
    ERR619514

Download the above sra files:

::

    mkdir sra_runs
    get-run.py -s -f sra_run_ids.txt -d sra_files

And export them to FASTQ files:

::

    mkdir fastqs
    for srr in $(cat sra_run_ids.txt)
    do
        fastq-dump                              \
            --split-spot                        \
            --skip-technical                    \
            --stdout                            \
            --readids                           \
            --defline-seq '@$sn/$ri'            \
            --defline-qual '+'                  \
            sra_files/${srr}.sra |               \
            gzip  > fastqs/${srr}.fastq.gz
    done

Don't worry about the full details of this ``fastq-dump``, but it produces a
gzipped interleaved fastq file.


Hashing
^^^^^^^

We directly utilise ``khmer``'s ``load-into-counting.py`` to hash reads to a
hash (Countgraph).

::

    mkdir hashes
    for srr in $(cat sra_run_ids.txt)
    do
        load-into-counting.py       \
            -N 1                    \
            -x 1e9                  \
            -k 20                   \
            -b                      \
            -f                      \
            -s tsv                  \
            hashes/${srr}.ct.gz
            fastqs/${srr}.fastq.gz
    done

This creates a hash with a single table and a billion bins for each run. Hashes
are saved, with gzip compression, to the ``*.ct.gz`` files under ``./hashes``.
These hashes are the direct input to ``kwip``. Note that this hash is probably
a bit small for this dataset, but we will go ahead anyway so this works on most
modern laptops.


Distance Calculation
^^^^^^^^^^^^^^^^^^^^

So here's the core of the protocol: calculating the pairwise distances between
these samples, which are from the two major groups of rice, Indica and
Japonica.

::

    kwip \
        -t 2 \
        -k rice.kern \
        -d rice.dist \
        hashes/*.ct.gz


This should calculate the weighted distance matrix between these samples, using
two threads.

Now, we plot these results using the R script ``img.R``. This creates plots of
the distance and kernel matrices, as well as a cluster dendrogram and
multi-dimensional scaling plot.

::

    Rscript img.R rice

This should create ``rice.pdf``. Inspect, and you should see two large
groupings corresponding to the two rice families.
