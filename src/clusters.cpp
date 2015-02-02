#include "clusters.hpp"


namespace ifgt
{


namespace
{


    // These file-local constants are used for initializing KMterm,
    // a magical structure that parameterizies the clustering.
    // The values are taken (with exceptions as noted) from kmlsample.cpp
    // in the kmlocal source. I guess we could optimize these values, but
    // that's a rabbit hole I'm not peeking into right now.
    static const double KMterm_a = 100;
    static const double KMterm_b = 0;
    static const double KMterm_c = 0;
    static const double KMterm_d = 0;
    static const double KMterm_mcr = 0.10;
    static const double KMterm_mar = 0.10;
    static const double KMterm_mrs = 3;
    static const double KMterm_ipa = 0.50; // beer o'clock?
    static const double KMterm_trl = 10;
    static const double KMterm_trf = 0.95;


}


Clusters::Clusters(const arma::mat& X, int K)
    : m_data(X.n_cols, X.n_rows)
    , m_centers(K, m_data)
    , m_term(KMterm_a, KMterm_b, KMterm_c, KMterm_d,
            KMterm_mcr, KMterm_mar, KMterm_mrs, KMterm_ipa,
            KMterm_trl, KMterm_trf)
    , m_algorithm(m_centers, m_term)
    , m_indices(X.n_rows)
    , m_rx()
{
    KMpointArray pt_array = m_data.getPts();
    for (size_t i = 0; i < X.n_rows; ++i)
    {
        for (size_t j = 0; j < X.n_cols; ++j)
        {
            pt_array[i][j] = KMcoord(X(i, j));
        }
    }
    m_data.buildKcTree();
    m_centers = m_algorithm.execute();

    std::vector<double> sq_dist(X.n_rows);
    m_centers.getAssignments(m_indices.data(), sq_dist.data());
    m_rx = std::sqrt(*std::max_element(sq_dist.begin(), sq_dist.end()));
}


}
