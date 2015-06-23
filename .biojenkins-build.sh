
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
rm -rf build
mkdir build
pushd build
    cmake $wdir -DKHMER_ROOT=$wdir/khmer/
    make
    ctest --verbose
popd
