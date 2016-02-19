#!/usr/bin/env sh
# Builds armadillo

set -ex

version=6.500.5
home=$(pwd)

wget https://downloads.sourceforge.net/project/arma/armadillo-$(version).tar.gz
tar xzf armadillo-$(version).tar.gz
rm armadillo-$(version).tar.gz
mv armadillo-$(version) armadillo
mkdir armadillo/build
cd armadillo/build
cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$(home)/local
make
make install
