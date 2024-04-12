#!/bin/bash

pwd
where cl.exe
export CC=cl.exe
export CXX=cl.exe
cmake .. \
    -G Ninja \
    -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
    -DWITH_TESTS=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -Dgtest_force_shared_crt=ON \
    -DBUILD_SHARED_LIBS=ON
