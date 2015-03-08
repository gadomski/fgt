#include "clustering.hpp"

#include "choose_truncation_number.hpp"
#include "constant_series.hpp"
#include "monomials.hpp"
#include "nchoosek.hpp"
#include "p_max_total.hpp"


namespace ifgt {


Clustering::Clustering(const arma::mat& source, int K, double bandwidth,
                       double epsilon)
    : m_source(source),
      m_indices(arma::zeros<arma::uvec>(source.n_rows)),
      m_centers(arma::zeros<arma::mat>(K, source.n_cols)),
      m_num_points(arma::zeros<arma::uvec>(K)),
      m_radii(K),
      m_rx(),
      m_bandwidth(bandwidth),
      m_epsilon(epsilon),
      m_p_max(0),
      m_constant_series(),
      m_is_initialized(false) {}


void Clustering::initialize() {
    m_p_max =
        choose_truncation_number(m_source.n_cols, m_bandwidth, m_epsilon, m_rx);
    m_constant_series = compute_constant_series(m_source.n_cols, m_p_max);
    m_is_initialized = true;
}


arma::mat Clustering::compute_C(const arma::vec& q) const {
    arma::mat C = arma::zeros<arma::mat>(
        m_centers.n_rows, get_p_max_total(m_source.n_cols, m_p_max));
    double h2 = m_bandwidth * m_bandwidth;

    for (arma::uword i = 0; i < m_source.n_rows; ++i) {
        arma::uword k = m_indices(i);
        arma::rowvec dx = m_source.row(i) - m_centers.row(k);
        double distance2 = arma::accu(arma::pow(dx, 2));
        arma::rowvec center_monomials =
            compute_monomials(dx / m_bandwidth, m_p_max);
        double f = q(i) * std::exp(-distance2 / h2);
        C.row(k) += f * center_monomials;
    }

    arma::rowvec constant_series =
        compute_constant_series(m_source.n_cols, m_p_max);
    for (arma::uword i = 0; i < C.n_rows; ++i) {
        C.row(i) = C.row(i) % constant_series;
    }

    return C;
}
}
