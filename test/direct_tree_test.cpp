#include "gtest/gtest.h"

#include "test/support.hpp"

namespace fgt {

TEST(DirectTree, MatchesDirect) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    auto expected = direct(source, target, 0.5);
    auto actual = direct_tree(source, target, 0.5, 1e-4);
    EXPECT_TRUE(actual.isApprox(expected, 1e-4));
}

TEST(DirectTree, WithWeights) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    Vector weights = Vector::LinSpaced(source.rows(), 0.1, 0.9);
    auto expected = direct(source, target, 0.5, weights);
    auto actual = direct_tree(source, target, 0.5, 1e-4, weights);
    EXPECT_TRUE(expected.isApprox(actual));
}

TEST(DirectTree, ClassBased) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    auto expected = direct(source, target, 0.5);
    auto actual = DirectTree(source, 0.5, 1e-4).compute(target);
    EXPECT_TRUE(actual.isApprox(expected, 1e-4));
}
}
