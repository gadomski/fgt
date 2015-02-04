#include <ifgt/clustering.hpp>

#include <iostream>

#include "choose_truncation_number.hpp"
#include "constant_series.hpp"
#include "monomials.hpp"
#include "nchoosek.hpp"
#include "p_max_total.hpp"


namespace ifgt
{


Clustering::Clustering(const arma::mat& X, int K, double h, double epsilon)
    : m_X(X)
    , m_indices(arma::zeros<arma::uvec>(X.n_rows))
    , m_centers(arma::zeros<arma::mat>(K, X.n_cols))
    , m_num_points(arma::zeros<arma::uvec>(K))
    , m_radii(K)
    , m_rx()
    , m_h(h)
    , m_epsilon(epsilon)
    , m_p_max(0)
    , m_constant_series()
{}


Clustering::~Clustering()
{}


void Clustering::compute()
{
    cluster();
    // TODO could check somehow to ensure a given clustring
    // populated everything
    m_p_max = choose_truncation_number(m_X.n_cols, m_h, m_epsilon, m_rx);
    m_constant_series = compute_constant_series(m_X.n_cols, m_p_max);
}


arma::mat Clustering::find_C(const arma::vec& q)
{
    arma::mat C = arma::zeros<arma::mat>(m_centers.n_rows, get_p_max_total(m_X.n_cols, m_p_max));
    double h2 = m_h * m_h;

    for (arma::uword i = 0; i < m_X.n_rows; ++i)
    {
        arma::uword k = m_indices(i);
        arma::rowvec dx = m_X.row(i) - m_centers.row(k);
        double distance2 = arma::accu(arma::pow(dx, 2));
        arma::rowvec center_monomials = compute_monomials(dx / m_h, m_p_max);
        double f = q(i) * std::exp(-distance2 / h2);
        C.row(k) += f * center_monomials;
    }

    arma::rowvec constant_series = compute_constant_series(m_X.n_cols, m_p_max);
    for (arma::uword i = 0; i < C.n_rows; ++i)
    {
        C.row(i) = C.row(i) % constant_series;
    }

    return C;
}


}
