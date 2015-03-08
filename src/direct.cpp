#include <fgt/fgt.hpp>


namespace fgt {


arma::vec Direct::compute_impl(const arma::mat& target,
                               const arma::vec& weights) const {
    double bandwidth2 = get_bandwidth() * get_bandwidth();
    const arma::mat& source = get_source();
    arma::vec g = arma::zeros<arma::vec>(target.n_rows);
    for (arma::uword j = 0; j < target.n_rows; ++j) {
        g(j) = arma::as_scalar(
            weights.t() *
            arma::exp(
                -arma::sum(arma::pow(source - arma::repmat(target.row(j),
                                                           source.n_rows, 1),
                                     2),
                           1) /
                bandwidth2));
    }
    return g;
}
}
