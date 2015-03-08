#include <ifgt/ifgt.hpp>

#include "clustering.hpp"
#include "clustering_factory.hpp"
#include "monomials.hpp"
#include "parameters.hpp"


namespace ifgt {


arma::vec ifgt(const arma::mat& source, const arma::mat& target,
               double bandwidth, double epsilon) {
    arma::vec q = arma::ones<arma::vec>(source.n_rows);
    return ifgt(source, target, bandwidth, epsilon, q);
}


arma::vec ifgt(const arma::mat& source, const arma::mat& target,
               double bandwidth, double epsilon, const arma::vec& q) {
    Parameters params = choose_parameters(source.n_cols, bandwidth, epsilon);
    return ifgt(source, target, bandwidth, epsilon, q, params);
}


arma::vec ifgt(const arma::mat& source, const arma::mat& target,
               double bandwidth, double epsilon, const arma::vec& q,
               const Parameters& params) {
    ClusteringFactory factory;
    Clustering clustering =
        factory.compute(source, params.K, bandwidth, epsilon);
    return ifgt(clustering, target, bandwidth, q, params);
}


arma::vec ifgt(const Clustering& clustering, const arma::mat& target,
               double bandwidth, const arma::vec& q, const Parameters& params) {
    // TODO check source.n_cols == target.n_cols
    arma::vec G(target.n_rows);
    arma::vec ry2 = arma::pow(params.r + clustering.get_radii(), 2);
    double h2 = bandwidth * bandwidth;
    arma::mat C = clustering.compute_C(q);
    for (arma::uword j = 0; j < target.n_rows; ++j) {
        G(j) = 0.0;
        for (arma::uword k = 0; k < params.K; ++k) {
            arma::rowvec dy = target.row(j) - clustering.get_centers().row(k);
            double distance2 = arma::accu(arma::pow(dy, 2));
            if (distance2 <= ry2(k)) {
                double g = std::exp(-distance2 / h2);
                G(j) += arma::accu(
                    C.row(k) %
                    compute_monomials(dy / bandwidth, clustering.get_p_max()) *
                    g);
            }
        }
    }
    return G;
}
}
