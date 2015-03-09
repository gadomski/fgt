#include "monomials.hpp"

#include "p_max_total.hpp"


namespace fgt {


void compute_monomials(arma::rowvec dx, arma::uword p_max,
                       std::vector<double>& monomials) {
    arma::vec heads = arma::zeros<arma::vec>(dx.n_cols);
    monomials[0] = 1;

    for (arma::uword k = 1, t = 1, tail = 1; k < p_max; ++k, tail = t) {
        for (arma::uword i = 0; i < dx.n_cols; ++i) {
            arma::uword head = heads[i];
            heads[i] = t;
            for (arma::uword j = head; j < tail; ++j, ++t) {
                monomials[t] = dx[i] * monomials[j];
            }
        }
    }
}
}
