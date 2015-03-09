#pragma once

#include <armadillo>


namespace fgt {


void compute_constant_series(arma::uword d, arma::uword p_max,
                             std::vector<double>& series);
}
