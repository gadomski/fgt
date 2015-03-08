#include "choose_truncation_number.hpp"

#include <gtest/gtest.h>


namespace ifgt
{


TEST(ChooseTruncationNumber, ReferenceImplementation)
{
    int d = 2;
    double h = 0.3;
    double epsilon = 1e-6;
    double rx = 0.1;
    int p_max = choose_truncation_number(d, h, epsilon, rx);
    EXPECT_EQ(9, p_max);
}
}
