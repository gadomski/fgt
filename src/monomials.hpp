#pragma once

#include <armadillo>


namespace fgt {


arma::rowvec compute_monomials(arma::rowvec dx, arma::uword p_max);
}
