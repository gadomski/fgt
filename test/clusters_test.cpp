#include <gtest/gtest.h>
#include "clusters.hpp"

#include <armadillo>

#include "config.hpp"


namespace ifgt
{


TEST(Clusters, Interface)
{
    arma::mat X;
    X.load(test_data_path("cluster-test.csv"));
    int K = 20;
    Clusters clusters = compute_clusters(X, K, 2);

    std::vector<arma::uword> expected_num_points =
        {287, 124, 107, 123, 514, 195, 122, 141, 410,
            288, 310, 406, 344, 402, 71, 273, 275, 254, 230, 124};
    std::vector<arma::uword> actual_num_points = arma::conv_to<std::vector<arma::uword>>::from(clusters.get_num_points());

    EXPECT_EQ(5000, clusters.get_indices().size());
    EXPECT_EQ(2, clusters.get_centers().n_cols);
    EXPECT_EQ(20, clusters.get_centers().n_rows);
    EXPECT_EQ(expected_num_points, actual_num_points);
    EXPECT_EQ(20, clusters.get_radii().size());
    EXPECT_NEAR(0.197607, clusters.get_rx(), 0.000001);
}


}
