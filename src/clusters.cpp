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



Clusters cluster(const arma::mat& X, int K)
{
    // I initially created a new KMdata constructor that accepted the
    // memptr() from the armadillo matrix. However, since kmlocal expects
    // row-major, but armadillo stores in colmn major, we were going to
    // have to copy the points *anyways*. So, we copy them once, explicity,
    // using the existing KMdata interface.
    KMdata km_pts(X.n_cols, X.n_rows);
    KMpointArray km_pt_array = km_pts.getPts();
    for (size_t i = 0; i < X.n_rows; ++i)
    {
        for (size_t j = 0; j < X.n_cols; ++j)
        {
            km_pts[i][j] = KMcoord(X(i, j));
        }
    }
    km_pts.buildKcTree();

    Clusters clusters = {
        KMfilterCenters(K, km_pts)
    };
    KMterm term(KMterm_a, KMterm_b, KMterm_c, KMterm_d, KMterm_mcr,
            KMterm_mar, KMterm_mrs, KMterm_ipa, KMterm_trl, KMterm_trf);
    // TODO do we want to allow the user to pick the algorithm? I don't see
    // why not.
    KMlocalHybrid algorithm(clusters.centers, term);

    clusters.centers = algorithm.execute();

    return clusters;
}


}
