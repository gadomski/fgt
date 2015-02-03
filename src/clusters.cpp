#include "clusters.hpp"

#include <iostream>


namespace ifgt
{


namespace
{


inline double ddist(const arma::rowvec& x, const arma::rowvec& y)
{
    return arma::sum(arma::pow(x - y, 2));
}


}


Clusters::Clusters(const arma::mat& X, int K)
    : m_X(X)
    , m_indices(X.n_rows)
    , m_centers(arma::zeros<arma::mat>(K, X.n_cols))
    , m_num_points(arma::zeros<arma::uvec>(K))
    , m_radii(K)
    , m_rx()
{}


void Clusters::compute()
{
    arma::uword N = m_indices.n_rows;
    arma::uword K = m_centers.n_rows;
    arma::uvec centers(K);
    arma::uvec cprev(N);
    arma::uvec cnext(N);
    arma::uvec far2c(K);
    arma::vec dist(N);

    std::random_device rd;
    //arma::uword nc = rd() % N;
    arma::uword nc = 2;
    centers[0] = nc;

    for (arma::uword i = 0; i < N; ++i)
    {
        dist(i) = (i == nc) ? 0.0 : ddist(m_X.row(i), m_X.row(nc));
        cnext[i] = i + 1;
        cprev[i] = i - 1;
    }
    
    cnext[N - 1] = 0;
    cprev[0] = N - 1;

    nc = std::max_element(dist.begin(), dist.end()) - dist.begin();
    far2c[0] = nc;
    m_radii[0] = dist(nc);

    for (int i = 1; i < K; ++i)
    {
        nc = far2c[std::max_element(m_radii.begin(), m_radii.begin() + i) - m_radii.begin()];

        centers(i) = nc;
        m_radii(i) = 0.0;
        dist(nc) = 0.0;
        m_indices(nc) = i;
        far2c(i) = nc;

        cnext(cprev(nc)) = cnext(nc);
        cprev(cnext(nc)) = cprev(nc);
        cnext(nc) = nc;
        cprev(nc) = nc;

        for (int j = 0; j < i; ++j)
        {
            arma::uword ct_j = centers(j);
            double dc2cq = ddist(m_X.row(ct_j), m_X.row(nc)) / 4;
            if (dc2cq < m_radii[j])
            {
                m_radii(j) = 0.0;
                far2c(j) = ct_j;
                arma::uword k = cnext(ct_j);
                while (k != ct_j)
                {
                    arma::uword nextk = cnext(k);
                    double dist2c_k = dist(k);
                    if (dc2cq < dist2c_k)
                    {
                        double dd = ddist(m_X.row(k), m_X.row(nc));
                        if (dd < dist2c_k)
                        {
                            dist(k) = dd;
                            m_indices(k) = i;
                            if (m_radii(i) < dd)
                            {
                                m_radii(i) = dd;
                                far2c(i) = k;
                            }
                            cnext(cprev(k)) = nextk;
                            cprev(nextk) = cprev(k);
                            cnext(k) = cnext(nc);
                            cprev(cnext(nc)) = k;
                            cnext(nc) = k;
                            cprev(k) = nc;
                        }
                        else if (m_radii(j) < dist2c_k)
                        {
                            m_radii(j) = dist2c_k;
                            far2c(j) = k;
                        }
                    }
                    else if (m_radii(j) < dist2c_k)
                    {
                        m_radii(j) = dist2c_k;
                        far2c(j) = k;
                    }
                    k = nextk;
                }
            }
        }
    }

    m_radii = arma::sqrt(m_radii);
    m_rx = *std::max_element(m_radii.begin(), m_radii.end());

    for (arma::uword i = 0; i < N; ++i)
    {
        ++m_num_points(m_indices(i));
        m_centers.row(m_indices(i)) += m_X.row(i);
    }

    m_centers = m_centers / arma::repmat(m_num_points, 1, 2);
}


Clusters compute_clusters(const arma::mat& X, int K)
{
    Clusters clusters(X, K);
    clusters.compute();
    return clusters;
}


}
