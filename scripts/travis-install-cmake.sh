#!/usr/bin/env sh
# Installs a pre-built cmake

set -ex

home=$(pwd)

wget https://cmake.org/files/v3.0/cmake-3.0.2.tar.gz
tar xzf cmake-3.0.2.tar.gz
rm cmake-3.0.2.tar.gz
cd cmake-3.0.2
./bootstrap
make
