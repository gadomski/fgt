#pragma once

#include <armadillo>


namespace ifgt {

class Clustering;

Clustering gonzalez_clustering(const arma::mat& source, int K, double bandwidth,
                               double epsilon, bool use_starting_idx,
                               arma::uword starting_idx);
}
