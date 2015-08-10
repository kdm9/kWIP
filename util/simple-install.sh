#!/bin/bash -xe

rm -rf build
mkdir -p build
cd build

cmake -DCMAKE_INSTALL_PREFIX=$HOME -DCMAKE_BUILD_TYPE=Release ..
make
make test

echo Everything worked!
echo Please cd into './build' and type make install
