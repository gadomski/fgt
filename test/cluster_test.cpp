#include "gtest/gtest.h"

#include "test/support.hpp"
#include "cluster.hpp"

namespace fgt {

TEST(Cluster, Reference) {
    size_t nclusters = 5;
    auto points =
        load_binary_test_matrix<double>("kmeans-128.dat", 128 * 128, 24);
    std::vector<double> starting_clusters(nclusters * points.cols);
    size_t index = 0;
    for (size_t j = 0; j < nclusters; ++j) {
        for (size_t k = 0; k < points.cols; ++k) {
            starting_clusters[j * points.cols + k] =
                points.data[index * points.cols + k];
        }
        index += 128 / nclusters;
    }
    Clustering clustering =
        cluster(points.data.data(), points.rows, points.cols, nclusters, 1e-5,
                starting_clusters);
    auto labels =
        load_binary_test_matrix<int>("kmeans-128-labels.dat", 128 * 128, 1);
    ASSERT_EQ(labels.data.size(), clustering.indices.size());
    for (size_t i = 0; i < labels.data.size(); ++i) {
        ASSERT_EQ(labels.data[i], clustering.indices[i]);
    }
}
}
