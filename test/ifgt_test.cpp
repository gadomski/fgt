#include "gtest/gtest.h"

#include "ifgt.hpp"
#include "test/support.hpp"

namespace fgt {

TEST(Ifgt, Reference) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    auto expected = direct(source, target, 0.5);
    auto actual = ifgt(source, target, 0.5, 1e-4);
    ASSERT_EQ(expected.size(), actual.size());
    EXPECT_TRUE(expected.isApprox(actual, 1e-2));
}

TEST(Ifgt, ClassBased) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    auto expected = direct(source, target, 0.5);
    auto actual = Ifgt(source, 0.5, 1e-4).compute(target);
    ASSERT_EQ(expected.size(), actual.size());
    EXPECT_TRUE(expected.isApprox(actual, 1e-2));
}

TEST(Ifgt, ChooseParameters) {
    IfgtParameters params = ifgt_choose_parameters(2, 0.3, 1e-6, 189, 200);
    EXPECT_EQ(13, params.nclusters);
    EXPECT_NEAR(1.1151, params.cutoff_radius, 1e-4);
}

TEST(Ifgt, ChooseTruncationNumber) {
    size_t truncation_number =
        ifgt_choose_truncation_number(2, 0.3, 1e-6, 0.1, 200);
    EXPECT_EQ(9, truncation_number);
}

TEST(Ifgt, HighBandwidth) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    target = target.array() + 2;
    auto expected = direct(source, target, 3.5);
    auto actual = ifgt(source, target, 3.5, 1e-4);
    ASSERT_EQ(expected.size(), actual.size());
    EXPECT_TRUE(expected.isApprox(actual, 1e-2));
}

TEST(Ifgt, UTM) {
    auto source = load_ascii_test_matrix("utm.txt");
    auto target = source;
    ASSERT_THROW(ifgt(source, target, 100, 1e-4), ifgt_no_clusters);
}

TEST(Ifgt, ManyDimensionsManyPoints) {
    Matrix source = Matrix::Random(10, 60);
    ASSERT_THROW(Ifgt(source, 0.4, 1e-4), fgt_error);
}
}
