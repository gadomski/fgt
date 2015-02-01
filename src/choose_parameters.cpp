#include "choose_parameters.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

#include "nchoosek.hpp"


namespace ifgt
{


namespace
{


    static const int TruncationNumberUpperLimit = 200;


}


Parameters choose_parameters(int d, double h, double epsilon, int k_limit)
{
    Parameters params;
    double R = std::sqrt(d);
    double h2 = h * h;
    double complexity_min = std::numeric_limits<double>::max();
    double rx = 0.0;

    params.r = std::min(R, h * std::sqrt(std::log(1 / epsilon)));
    params.k = 1;
    params.p_max = 0;

    for (int i = 0; i < k_limit; ++i)
    {
        rx = std::pow(double(i + 1), -1.0 / double(d));
        double rx2 = rx * rx;
        double n = std::min(double(i + 1), std::pow(params.r / rx, double(d)));
        double error = 1;
        double temp = 1;
        int p = 0;

        while ((error > epsilon) and (p <= TruncationNumberUpperLimit))
        {
            ++p;
            double b = std::min((rx + std::sqrt(rx2 + 2 * p * h2)) / 2,
                rx + params.r);
            double c = rx - b;
            temp *= 2 * rx * b / h2 / p;
            error = temp * std::exp(- c * c / h2);
        }
        double complexity = i + 1 + std::log(double(i + 1)) + (i + 1) * nchoosek(p - 1 + d, d);

        if (complexity < complexity_min)
        {
            complexity_min = complexity;
            params.k = i + 1;
            params.p_max = p;
        }
    }

    return params;
}


}
