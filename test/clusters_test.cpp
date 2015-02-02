#include <gtest/gtest.h>
#include "clusters.hpp"

#include <armadillo>


namespace ifgt
{


TEST(Clusters, Interface)
{
    arma::mat X = arma::randu<arma::mat>(5000, 2);
    int K = 20;
    Clusters clusters = cluster(X, K);
    // TODO test out some of the characteristics of a cluster
}


}
