#include "clustering/gonzalez.hpp"

#include "clustering_factory.hpp"
#include "clustering.hpp"
#include "config.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace ifgt
{


TEST(GonzalezClustering, ReferenceImplementation)
{
    arma::mat X;
    X.load(test_data_path("X.csv"));
    int K = 15;
    double h = 0.4;
    double epsilon = 1e-3;

    Clustering clustering = gonzalez_clustering(X, K, h, epsilon, true, 2);
    clustering.initialize();

    std::vector<arma::uword> expected_num_points = {167, 167, 185, 470, 482,
                                                    168, 177, 179, 571, 356,
                                                    825, 167, 417, 362, 307};
    std::vector<arma::uword> actual_num_points =
        arma::conv_to<std::vector<arma::uword>>::from(
            clustering.get_num_points());

    EXPECT_EQ(5000, clustering.get_indices().size());
    EXPECT_EQ(2, clustering.get_centers().n_cols);
    EXPECT_EQ(15, clustering.get_centers().n_rows);
    EXPECT_EQ(expected_num_points, actual_num_points);
    EXPECT_EQ(15, clustering.get_radii().size());
    EXPECT_NEAR(0.0799, clustering.get_radius(0), 0.0001);
    EXPECT_NEAR(0.1838, clustering.get_rx(), 0.0001);
    EXPECT_TRUE(clustering.get_p_max() > 0);
}
}
