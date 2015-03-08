#include <fgt/fgt.hpp>

#include "clustering.hpp"
#include "clustering/gonzalez.hpp"
#include "monomials.hpp"
#include "parameters.hpp"


namespace fgt {


Ifgt::Ifgt(const arma::mat& source, double bandwidth, double epsilon,
           int k_limit)
    : GaussianTransform(source, bandwidth),
      m_epsilon(epsilon),
      m_k_limit(k_limit),
      m_clustering_starting_index(std::make_pair(false, 0)) {}


optional_arma_uword_t Ifgt::get_clustering_starting_index() const {
    return m_clustering_starting_index;
}


Ifgt& Ifgt::set_clustering_starting_index(arma::uword index) {
    m_clustering_starting_index = std::make_pair(true, index);
    return *this;
}


arma::vec Ifgt::compute(const arma::mat& target,
                        const arma::vec& weights) const {
    const arma::mat& source = get_source();
    double bandwidth = get_bandwidth();
    Parameters params = choose_parameters(source.n_cols, bandwidth, m_epsilon);
    Clustering clustering =
        gonzalez_clustering(source, params.K, bandwidth, m_epsilon,
                            get_clustering_starting_index());
    // TODO check source.n_cols == target.n_cols
    arma::vec G(target.n_rows);
    arma::vec ry2 = arma::pow(params.r + clustering.get_radii(), 2);
    double h2 = bandwidth * bandwidth;
    arma::mat C = clustering.compute_C(weights);
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
