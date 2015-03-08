#include <fgt/fgt.hpp>

#include "config.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt {


TEST(Direct, ReferenceImplementation) {
    arma::mat source, target;
    source.load(test_data_path("X.csv"));
    target.load(test_data_path("Y.csv"));
    source = source.rows(0, 99);
    target = target.rows(0, 99);
    double bandwidth = 0.4;

    Direct direct(source, bandwidth);
    arma::vec g = direct.compute(target);

    EXPECT_EQ(g.n_rows, 100);
    EXPECT_DOUBLE_EQ(4.3701728529168333, arma::min(g));
    EXPECT_DOUBLE_EQ(93.72376652714631, arma::max(g));
}
}
