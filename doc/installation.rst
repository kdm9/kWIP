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
releases page <https://github.com/kdmurray91/kWIP/releases>`_ in an archive
named ``kwip-binaries_VERSION.tar.xz`` or similar (where ``VERSION`` is the
latest released version). Please download this to your **GNU/Linux** machine,
and:

.. code-block:: shell

    # Assuming you want to install to ~/bin/kwip
    PREFIX=$HOME
    cd $PREFIX

    # If ~/bin/ is not in $PATH, you won't be able to use kwip.
    # Perform the command below to ensure PATH is set correctly.
    echo "PATH=\"${PREFIX}/bin:\${PATH}\"" >> ~/.bashrc
    . ~/.bashrc

    # Below, repace all text in quotes with the URL to the archive
    # on the GitHub release page
    wget "URL FROM THE RELEASES PAGE"

    # Extract the precompile binaries and documentation
    tar xvf kwip-binaries*.tar.xz

    # Check your installation by typing
    kwip --help

    # The HTML documenation is available under ./share/doc/kwip
    ls $PREFIX/share/doc/kwip
