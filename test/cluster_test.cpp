#include "gtest/gtest.h"

#include "test/support.hpp"
#include "cluster.hpp"

namespace fgt {

TEST(Cluster, Reference) {
    Matrix::Index nclusters = 5;
    auto points =
        load_binary_test_matrix<double>("kmeans-128.dat", 128 * 128, 24);
    Matrix starting_clusters(nclusters, points.cols());
    Matrix::Index index = 0;
    auto N = points.cols();
    for (Matrix::Index j = 0; j < nclusters; ++j) {
        for (Matrix::Index k = 0; k < N; ++k) {
            starting_clusters(j, k) = points(index, k);
        }
        index += 128 / nclusters;
    }
    Clustering clustering = cluster(points, nclusters, 1e-5, starting_clusters);
    auto labels =
        load_binary_test_matrix<int>("kmeans-128-labels.dat", 128 * 128, 1);
    ASSERT_EQ(labels.size(), clustering.indices.size());
    EXPECT_TRUE(labels.cast<Matrix::Index>().isApprox(clustering.indices));
}
}
