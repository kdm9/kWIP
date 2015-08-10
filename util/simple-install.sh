#!/bin/bash -xe

rm -rf build
mkdir -p build

cd build

cmake -DCMAKE_INSTALL_PREFIX=$HOME -DCMAKE_BUILD_TYPE=Release ..

echo OK, now cd into './build' and type make and/or make install
