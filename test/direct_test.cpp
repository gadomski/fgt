#include <fgt/direct.hpp>

#include "config.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt {


TEST(Direct, ReferenceImplementation) {
    arma::mat X, Y;
    X.load(test_data_path("X.csv"));
    Y.load(test_data_path("Y.csv"));
    X = X.rows(0, 99);
    Y = Y.rows(0, 99);
    double h = 0.4;

    arma::vec G = direct(X, Y, h);

    EXPECT_EQ(G.n_rows, 100);
    EXPECT_DOUBLE_EQ(4.3701728529168333, arma::min(G));
    EXPECT_DOUBLE_EQ(93.72376652714631, arma::max(G));
}
}
