#include "gtest/gtest.h"

#include "ifgt.hpp"
#include "test/support.hpp"

namespace fgt {

TEST(Ifgt, Reference) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    double bandwidth = 0.5;
    double epsilon = 1e-4;
    auto expected = direct(source, target, bandwidth);
    auto actual = ifgt(source, target, bandwidth, epsilon);
    ASSERT_EQ(expected.size(), actual.size());
    EXPECT_LT((expected - actual).array().abs().maxCoeff() / actual.size(), epsilon);
}

TEST(Ifgt, ClassBased) {
    auto source = load_ascii_test_matrix("X.txt");
    auto target = load_ascii_test_matrix("Y.txt");
    double bandwidth = 0.5;
    double epsilon = 1e-4;
    auto expected = direct(source, target, bandwidth);
    auto actual = Ifgt(source, bandwidth, epsilon).compute(target);
    ASSERT_EQ(expected.size(), actual.size());
    EXPECT_LT((expected - actual).array().abs().maxCoeff() / actual.size(), epsilon);
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
    double bandwidth = 3.5;
    double epsilon = 1e-4;
    auto expected = direct(source, target, bandwidth);
    auto actual = ifgt(source, target, bandwidth, epsilon);
    ASSERT_EQ(expected.size(), actual.size());
    EXPECT_LT((expected - actual).array().abs().maxCoeff() / actual.size(), epsilon);
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
