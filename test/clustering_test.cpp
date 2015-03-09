#include "clustering.hpp"

#include "config.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt {


TEST(GonzalezClustering, ReferenceImplementation) {
    arma::mat X;
    X.load(test_data_path("X.csv"));
    int K = 15;
    double h = 0.4;
    double epsilon = 1e-3;

    GonzalezClustering clustering(X, K, h, epsilon, {true, 2});
    clustering.cluster();

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
    EXPECT_DOUBLE_EQ(0.079945937357691899, clustering.get_radius(0));
    EXPECT_DOUBLE_EQ(0.18376119993077969, clustering.get_rx());
    EXPECT_EQ(6, clustering.get_p_max());
}
}
