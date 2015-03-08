#include <fgt/direct.hpp>

#include "config.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt {


TEST(Direct, ReferenceImplementation) {
    arma::mat X, Y;
    X.load(test_data_path("X.csv"));
    Y.load(test_data_path("Y.csv"));
    double h = 0.4;

    arma::vec G = direct(X, Y, h);

    EXPECT_EQ(G.n_rows, 5000);
    EXPECT_DOUBLE_EQ(359.47775885562908, arma::min(G));
    EXPECT_DOUBLE_EQ(2221.3622088016318, arma::max(G));
}
}
