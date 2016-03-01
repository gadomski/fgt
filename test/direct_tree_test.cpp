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
    Vector weights = Vector::Ones(source.rows());
    auto no_weights = direct_tree(source, target, 0.5, 1e-4);
    auto with_weights = direct_tree(source, target, 0.5, 1e-4, weights);
    EXPECT_TRUE(no_weights.isApprox(with_weights));
}

TEST(DirectTree, ClassBased) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    auto expected = direct(source, target, 0.5);
    auto actual = DirectTree(source, 0.5, 1e-4).compute(target);
    EXPECT_TRUE(actual.isApprox(expected, 1e-4));
}
}
