#include <fgt/fgt.hpp>

#include "config.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt {


TEST(Direct, OneTargetOneSource) {
    arma::mat source = {0, 0};
    source.reshape(1, 2);
    arma::mat target = {1, 2};
    target.reshape(1, 2);

    Direct direct(source, 1);
    arma::vec g = direct.compute(target);
    EXPECT_EQ(g.n_rows, 1);
    EXPECT_DOUBLE_EQ(0.006737946999085467, g(0));
}


TEST(Direct, TwoTargetsTwoSources) {
    arma::mat source = {0, 0, 0, 1};
    source.reshape(2, 2);
    arma::mat target = {0, 0.5, 0.25, 0.25};
    target.reshape(2, 2);

    Direct direct(source, 0.2);
    arma::vec g = direct.compute(target);
    EXPECT_EQ(g.n_rows, 2);
    EXPECT_DOUBLE_EQ(0.2096121683000387, g(0));
    EXPECT_DOUBLE_EQ(0.00040464667729846903, g(1));
}
}
