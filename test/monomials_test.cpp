#include "monomials.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace ifgt
{


TEST(Monomials, ReferenceImplementation)
{
    arma::rowvec dx = {0.319529, 0.401860};
    arma::uword p_max = 70;
    arma::rowvec monomials = compute_monomials(dx / 0.4, p_max);
    EXPECT_EQ(2485, monomials.size());
    EXPECT_NEAR(1.0, monomials.at(0), 0.000001);
    EXPECT_NEAR(0.798823, monomials.at(1), 0.000001);
    EXPECT_NEAR(1.004649, monomials.at(2), 0.000001);
    EXPECT_NEAR(0.638118, monomials.at(3), 0.000001);
    EXPECT_NEAR(0.802537, monomials.at(4), 0.000001);
    EXPECT_NEAR(1.009321, monomials.at(5), 0.000001);
}
}
