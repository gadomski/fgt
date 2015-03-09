#include "truncation_number.hpp"

#include <cmath>


namespace fgt {


arma::uword choose_truncation_number(double distance2, double cutoff_radius,
                                     double bandwidth, double epsilon,
                                     arma::uword max_truncation_number) {
    double distance = std::sqrt(distance2);
    double bandwidth2 = bandwidth * bandwidth;
    double error = epsilon + 1;
    double temp = 1;
    arma::uword truncation_number = 0;

    while (error > epsilon and truncation_number < max_truncation_number) {
        ++truncation_number;
        double b = std::min(
            (distance +
             std::sqrt(distance2 + 2 * truncation_number * bandwidth2)) /
                2,
            cutoff_radius);
        double c = distance - b;
        temp *= 2 * distance * b / bandwidth2 / double(truncation_number);
        error = temp * std::exp(-c * c / bandwidth2);
    }
    return truncation_number;
}
}
