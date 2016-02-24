#include "gtest/gtest.h"

#include "test/support.hpp"

namespace fgt {

TEST(Direct, Reference) {
    auto source = load_ascii_test_matrix<double>("X.txt");
    auto target = load_ascii_test_matrix<double>("Y.txt");
    auto actual = direct(source.data.data(), source.rows, target.data.data(),
                         target.rows, source.cols, 0.5);
    auto expected = load_ascii_test_matrix<double>("direct.txt").data;
    for (size_t i = 0; i < actual.size(); ++i) {
        ASSERT_NEAR(expected[i], actual[i], 1e-4);
    }
}

TEST(Direct, WithWeights) {
    auto source = load_ascii_test_matrix<double>("X.txt");
    auto target = load_ascii_test_matrix<double>("Y.txt");
    auto no_weights = direct(source.data.data(), source.rows,
                             target.data.data(), target.rows, source.cols, 0.5);
    std::vector<double> weights(source.rows, 1.0);
    auto with_weights =
        direct(source.data.data(), source.rows, target.data.data(), target.rows,
               source.cols, 0.5, weights.data());
    for (size_t i = 0; i < no_weights.size(); ++i) {
        ASSERT_DOUBLE_EQ(no_weights[i], with_weights[i]);
    }
}

TEST(Direct, ClassBased) {
    auto source = load_ascii_test_matrix<double>("X.txt");
    auto target = load_ascii_test_matrix<double>("Y.txt");
    Direct direct(source.data.data(), source.rows, source.cols, 0.5);
    auto actual = direct.compute(target.data.data(), target.rows);
    auto expected = load_ascii_test_matrix<double>("direct.txt").data;
    for (size_t i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], actual[i], 1e-4);
    }
}
}
