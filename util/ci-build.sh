set -xe

export wdir=`pwd`

BUILD_TYPE="${BUILD_TYPE:-Debug}"

# Build kwip
rm -rf build target
mkdir build
mkdir target
pushd build
    cmake \
        $wdir \
        -DCMAKE_BUILD_TYPE=$BUILD_TYPE \
        -DCMAKE_INSTALL_PREFIX=$wdir/target
    make
    ctest --verbose
    make install
    test -x $wdir/target/bin/kwip
popd


# Build documentation

if [ -x "$(which sphinx-build)" ]
then
    pushd doc
        make
    popd
fi
