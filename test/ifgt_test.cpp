#include "gtest/gtest.h"

#include "test/support.hpp"
#include "ifgt.hpp"

namespace fgt {

TEST(Ifgt, Reference) {
    auto source = load_ascii_test_matrix<double>("X.txt");
    auto target = load_ascii_test_matrix<double>("Y.txt");
    auto expected = direct(source.data.data(), source.rows, target.data.data(),
                           target.rows, source.cols, 0.5);
    auto actual = ifgt(source.data.data(), source.rows, target.data.data(),
                       target.rows, source.cols, 0.5, 1e-4);
    ASSERT_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], actual[i], 1e-2);
    }
}

TEST(Ifgt, ClassBased) {
    auto source = load_ascii_test_matrix<double>("X.txt");
    auto target = load_ascii_test_matrix<double>("Y.txt");
    auto expected = direct(source.data.data(), source.rows, target.data.data(),
                           target.rows, source.cols, 0.5);
    auto actual = Ifgt(source.data.data(), source.rows, source.cols, 0.5, 1e-4)
                      .compute(target.data.data(), target.rows);
    ASSERT_EQ(expected.size(), actual.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        ASSERT_NEAR(expected[i], actual[i], 1e-2);
    }
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
}
