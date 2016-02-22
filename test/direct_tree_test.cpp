#include "gtest/gtest.h"

#include "test/support.hpp"

namespace fgt {

TEST(DirectTree, MatchesDirect) {
    auto source = load_ascii_test_matrix<double>("X.txt");
    auto target = load_ascii_test_matrix<double>("Y.txt");
    auto expected = direct(source.data.data(), source.rows, target.data.data(),
                           target.rows, source.cols, 0.5);
    auto actual =
        direct_tree(source.data.data(), source.rows, target.data.data(),
                    target.rows, source.cols, 0.5, 1e-4);
    for (size_t i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], actual[i], 1e-4);
    }
}

TEST(DirectTree, WithWeights) {
    auto source = load_ascii_test_matrix<double>("X.txt");
    auto target = load_ascii_test_matrix<double>("Y.txt");
    std::vector<double> weights(source.rows, 1.0);
    auto no_weights =
        direct_tree(source.data.data(), source.rows, target.data.data(),
                    target.rows, source.cols, 0.5, 1e-4);
    auto with_weights =
        direct_tree(source.data.data(), source.rows, target.data.data(),
                    target.rows, source.cols, 0.5, 1e-4, weights.data());
    for (size_t i = 0; i < no_weights.size(); ++i) {
        ASSERT_DOUBLE_EQ(no_weights[i], with_weights[i]);
    }
}
}
