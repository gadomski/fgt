#include "constant_series.hpp"

#include "p_max_total.hpp"


namespace fgt {


void compute_constant_series(arma::uword d, arma::uword p_max,
                                     std::vector<double>& series) {
    arma::uvec heads = arma::zeros<arma::uvec>(d + 1);
    heads(d) = std::numeric_limits<arma::uword>::max();
    arma::uvec cinds = arma::zeros<arma::uvec>(series.size());
    series[0] = 1;

    for (arma::uword k = 1, t = 1, tail = 1; k < p_max; ++k, tail = t) {
        for (arma::uword i = 0; i < heads.n_rows; ++i) {
            arma::uword head = heads(i);
            heads(i) = t;
            for (arma::uword j = head; j < tail; ++j, ++t) {
                cinds(t) = (j < heads(i + 1)) ? (cinds(j) + 1) : 1;
                series[t] = 2.0 * series[j] / double(cinds(t));
            }
        }
    }
}
}
