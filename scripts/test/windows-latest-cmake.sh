#!/bin/bash


export PATH=`pwd`/../bin;$PATH
ctest -VV --output-on-failure
