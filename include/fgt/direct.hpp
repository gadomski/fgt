#pragma once

#include <armadillo>


namespace fgt {


arma::vec direct(const arma::mat& source, const arma::mat& target,
                 double bandwidth);
arma::vec direct(const arma::mat& source, const arma::mat& target,
                 double bandwidth, const arma::vec& weights);
}
