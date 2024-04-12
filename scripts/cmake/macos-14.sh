#!/bin/bash

cmake .. \
    -G Ninja \
    -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX \
    -DWITH_TESTS=ON \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_SHARED_LIBS=ON
