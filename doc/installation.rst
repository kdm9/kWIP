===================
Installing ``kWIP``
===================

The full way
------------

Dependencies:

- ``zlib``
- ``cmake>=2.8``
- Optionally, ``liboxli`` and ``Eigen3`` are required. These libraries are bundled
  with kWIP, and the internal copy will be used if system copies are not.
- A C++11 compiler that supports OpenMP (i.e. gcc >=4.8)

On Debian (or Debian derivatives) the dependencies of ``kWIP`` can be installed
with:

.. code-block:: shell

    sudo apt-get install zlib1g-dev cmake build-essential git

Then, to compile ``kWIP``:

.. code-block:: shell

    git clone https://github.com/kdmurray91/kWIP.git
    cd kWIP
    mkdir build && cd build
    cmake .. -DCMAKE_INSTALL_PREFIX=${HOME}
    make
    make test
    make install

The commands above assume you want to install kWIP to your home directory. This
is probably required on clusters, and necessary without root privileges. To
install to, e.g, ``/usr/local/``, replace all occurrences of ``$HOME`` with your
preferred installation prefix.

The easy way
------------

Pre-compiled static binaries for 64-bit GNU/Linux are provided on the `GitHub
releases page <https://github.com/kdmurray91/kWIP/releases>`_. The following
commands will obtain and install ``kWIP``.

.. code-block:: shell

    # If ~/bin/ is not in $PATH, you won't be able to use kwip.
    # Perform the command below to ensure PATH is set correctly.
    echo "PATH=\"${HOME}/bin:\${PATH}\"" >> ~/.bashrc
    . ~/.bashrc

    cd $HOME
    wget https://github.com/kdmurray91/kWIP/releases/download/0.2.0/kwip-binaries_0.2.0.tar.gz
    wget https://github.com/kdmurray91/kWIP/releases/download/0.2.0/kwip-binaries_0.2.0.tar.gz.sha256sums
    sha256sum -c kwip-binaries_0.2.0.tar.gz.sha256sums

    # Extract the precompiled binaries
    tar xvf kwip-binaries*.tar.gz

    # Check your installation by typing
    kwip --help


Please note that as these binaries are compiled to be as widely compatible as
possible, the compiler will use few modern optimisations. Therefore it is
possible that these static binaries will be slower on modern processors than
binaries compiled from source.
