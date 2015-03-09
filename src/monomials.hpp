#pragma once

#include <armadillo>


namespace fgt {


void compute_monomials(arma::rowvec dx, arma::uword p_max,
                       std::vector<double>& monomials);
}
