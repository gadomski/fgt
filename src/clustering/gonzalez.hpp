#pragma once

#include <fgt/typedefs.hpp>

#include <armadillo>


namespace fgt {

class Clustering;

Clustering gonzalez_clustering(const arma::mat& source, int K, double bandwidth,
                               double epsilon, optional_arma_uword_t starting_index);
}
