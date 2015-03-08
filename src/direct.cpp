#include <fgt/direct.hpp>


namespace fgt {


arma::vec direct(const arma::mat& source, const arma::mat& target,
                 double bandwidth) {
    arma::vec weights = arma::ones<arma::vec>(source.n_rows);
    return direct(source, target, bandwidth, weights);
}


arma::vec direct(const arma::mat& source, const arma::mat& target,
                 double bandwidth, const arma::vec& weights) {
    double bandwidth2 = bandwidth * bandwidth;
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
