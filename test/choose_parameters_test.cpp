#include <gtest/gtest.h>
#include "choose_parameters.hpp"


namespace ifgt
{


TEST(ChooseParameters, ImplementationExample)
{
    int d = 2;
    double h = 0.3;
    double epsilon = 1e-6;
    int k_limit = 189;
    Parameters params = choose_parameters(d, h, epsilon, k_limit);
    EXPECT_EQ(13, params.k);
    EXPECT_EQ(17, params.p_max);
    EXPECT_NEAR(1.1151, params.r, 0.0001);
}


}
