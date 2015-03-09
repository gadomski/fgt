#include "truncation_number.hpp"

#include <gtest/gtest.h>

#include <cmath>


namespace fgt
{


TEST(ChooseTruncationNumber, ReferenceImplementation) {
    double bandwidth = 0.3;
    double epsilon = 1e-6;
    double distance2 = 0.01;
    double cutoff_radius = std::sqrt(distance2) + bandwidth * std::sqrt(std::log(1 / epsilon));
    arma::uword max_truncation_number = 300;
    int p_max = choose_truncation_number(distance2, cutoff_radius, bandwidth, epsilon);
    EXPECT_EQ(9, p_max);
}


}
