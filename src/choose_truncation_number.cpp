#include "choose_truncation_number.hpp"

#include <algorithm>
#include <cmath>


namespace ifgt
{


namespace
{
// TODO I don't know why this number differs from that in
// choose_parameters. We should probably consolidate.
static const int TruncationNumberUpperLimit = 300;
}


int choose_truncation_number(int d, double bandwidth, double epsilon, double rx)
{
    double r =
        std::min(std::sqrt(d), bandwidth * std::sqrt(std::log(1 / epsilon)));
    double rx2 = rx * rx;
    double h2 = bandwidth * bandwidth;
    double error = 1;
    double temp = 1;
    int p = 0;

    while ((error > epsilon) and (p <= TruncationNumberUpperLimit))
    {
        ++p;
        double b = std::min((rx + std::sqrt(rx2 + 2 * p * h2)) / 2, rx + r);
        double c = rx - b;
        temp *= 2 * rx * b / h2 / p;
        error = temp * std::exp(-c * c / h2);
    }

    return p;
}
}
