#include "constant_series.hpp"

#include "p_max_total.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt {


TEST(ConstantSeries, ReferenceImplementation) {
    arma::uword d = 2;
    arma::uword p_max = 70;
    std::vector<double> series(get_p_max_total(d, p_max));
    compute_constant_series(2, 70, series);
    EXPECT_EQ(2485, series.size());
    EXPECT_DOUBLE_EQ(1, series.at(0));
    EXPECT_DOUBLE_EQ(2, series.at(1));
    EXPECT_DOUBLE_EQ(2, series.at(2));
    EXPECT_DOUBLE_EQ(2, series.at(3));
    EXPECT_DOUBLE_EQ(4, series.at(4));
    EXPECT_DOUBLE_EQ(2, series.at(5));
    EXPECT_DOUBLE_EQ(4.0 / 3.0, series.at(6));
}
}
