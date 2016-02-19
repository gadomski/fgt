#!/usr/bin/env sh
# Runs a suite of benchmarks on the fgt C++ library.

modes="direct direct_tree ifgt"
rows="1000 5000 10000"
cols="3"
bandwidths="1 0.562 0.316 0.178 0.1"
niter=5

for mode in ${modes}; do
    for row in ${rows}; do
        for col in ${cols}; do
            for bandwidth in ${bandwidths}; do
                printf "%s %s %s %s" $mode $row $col $bandwidth
                for i in `seq 1 ${niter}`; do
                    printf " %s" $(build/bench/bench $mode random $row $col $bandwidth)
                done
                printf "\n"
            done
        done
    done
done
