#include "gtest/gtest.h"

#include "test/support.hpp"

namespace fgt {

TEST(DirectTree, MatchesDirect) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    double bandwidth = 0.5;
    auto expected = direct(source, target, bandwidth);
    double epsilon = 1e-4;
    auto actual = direct_tree(source, target, bandwidth, epsilon);
    EXPECT_LT((expected - actual).array().abs().maxCoeff() / actual.size(), epsilon);
}

TEST(DirectTree, WithWeights) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    Vector weights = Vector::LinSpaced(source.rows(), 0.1, 0.9);
    double bandwidth = 0.5;
    auto expected = direct(source, target, bandwidth, weights);
    double epsilon = 1e-4;
    auto actual = direct_tree(source, target, bandwidth, epsilon, weights);
    EXPECT_LT((expected - actual).array().abs().maxCoeff() / weights.sum(), epsilon);
}

TEST(DirectTree, ClassBased) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    double bandwidth = 0.5;
    auto expected = direct(source, target, bandwidth);
    double epsilon = 1e-4;
    auto actual = DirectTree(source, bandwidth, epsilon).compute(target);
    EXPECT_LT((expected - actual).array().abs().maxCoeff() / actual.size(), epsilon);
}

TEST(DirectTree, Fish) {
    auto source = load_ascii_test_matrix("fish.txt");
    auto target = load_ascii_test_matrix("fish-transformed.txt");
    double bandwidth = 1.0;
    auto expected = direct(source, target, bandwidth);
    auto actual = direct_tree(source, target, bandwidth, 1e-5);
    EXPECT_TRUE(expected.isApprox(actual, 1e-4));
}
}
