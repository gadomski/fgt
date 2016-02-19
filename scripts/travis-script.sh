#!/usr/bin/env sh
# Travis's build and test script

set -ex

home=$(pwd)

mkdir build
cd build
cmake .. -DWITH_TESTS=ON \
  -DWITH_OPENMP=${FGT_WITH_OPENMP} \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_VERBOSE_MAKEFILE=${FGT_CMAKE_VERBOSE_MAKEFILE} \
  -DCMAKE_PREFIX_PATH=${home}/local \
  -DCMAKE_INSTALL_PREFIX=${home}/local
make
make test
make install
