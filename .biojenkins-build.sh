
set -xe

export wdir=`pwd`

test -d khmer-src || git clone https://github.com/ged-lab/khmer khmer-src

# update and install khmer
pushd khmer-src
    git fetch --all
    git reset --hard origin/master
    rm -rf $wdir/khmer
    mkdir $wdir/khmer
    pushd lib
        make install PREFIX=$wdir/khmer/
    popd
popd


# Build kwip
rm -rf build target
mkdir build
mkdir target
pushd build
    cmake $wdir -DKHMER_ROOT=$wdir/khmer/ -DCMAKE_INSTALL_PREFIX=$wdir/target
    make
    ctest --verbose
    make install
    test -x $wdir/target/bin/kwip
popd
