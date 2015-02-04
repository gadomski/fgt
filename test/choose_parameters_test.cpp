#include <gtest/gtest.h>
#include <ifgt/choose_parameters.hpp>

#include <armadillo>


namespace ifgt
{


TEST(ChooseParameters, ImplementationExample)
{
    int d = 2;
    double h = 0.3;
    double epsilon = 1e-6;
    int k_limit = 189;
    Parameters params = choose_parameters(d, h, epsilon, k_limit);
    EXPECT_EQ(13, params.K);
    EXPECT_EQ(17, params.p_max);
    EXPECT_NEAR(1.1151, params.r, 0.0001);
}


TEST(ChooseParameters, ReferenceImplementation)
{
    arma::uword d = 2;
    double h = 0.4;
    double epsilon = 1e-3;
    arma::uword k_limit = 50;
    Parameters params = choose_parameters(d, h, epsilon, k_limit);
    EXPECT_EQ(15, params.K);
    EXPECT_EQ(8, params.p_max);
    EXPECT_NEAR(1.051304, params.r, 0.000001);
}


}
