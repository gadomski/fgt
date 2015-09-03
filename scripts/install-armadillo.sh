#!/usr/bin/env sh
ARMADILLO_VERSION=4.650.3
INSTALL_PREFIX=$1

set -e
wget https://downloads.sourceforge.net/project/arma/armadillo-$ARMADILLO_VERSION.tar.gz --no-check-certificate
tar -xvzf armadillo-$ARMADILLO_VERSION.tar.gz
mkdir armadillo-$ARMADILLO_VERSION/build
cd armadillo-$ARMADILLO_VERSION/build && \
    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$INSTALL_PREFIX && \
    make && \
    make install
