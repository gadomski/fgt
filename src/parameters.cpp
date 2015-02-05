#include <ifgt/parameters.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

#include <armadillo>

#include "nchoosek.hpp"


namespace ifgt
{


namespace
{


    static const arma::uword TruncationNumberUpperLimit = 200;
    static const arma::uword KLimitFactor = 20;

}


Parameters choose_parameters(arma::uword d, double bandwidth, double epsilon)
{
    return choose_parameters(d, bandwidth, epsilon, std::round(KLimitFactor / bandwidth));
}


Parameters choose_parameters(arma::uword d, double bandwidth, double epsilon, arma::uword k_limit)
{
    Parameters params;
    double R = std::sqrt(d);
    double h2 = bandwidth * bandwidth;
    double complexity_min = std::numeric_limits<double>::max();
    double rx = 0.0;

    params.r = std::min(R, bandwidth * std::sqrt(std::log(1 / epsilon)));
    params.K = 1;

    for (arma::uword i = 0; i < k_limit; ++i)
    {
        rx = std::pow(double(i + 1), -1.0 / double(d));
        double rx2 = rx * rx;
        double n = std::min(double(i + 1), std::pow(params.r / rx, double(d)));
        double error = 1;
        double temp = 1;
        arma::uword p = 0;

        while ((error > epsilon) and (p <= TruncationNumberUpperLimit))
        {
            ++p;
            double b = std::min((rx + std::sqrt(rx2 + 2 * p * h2)) / 2,
                rx + params.r);
            double c = rx - b;
            temp *= 2 * rx * b / h2 / p;
            error = temp * std::exp(- c * c / h2);
        }
        double complexity = i + 1 + std::log(double(i + 1)) + (n + 1) * nchoosek(p - 1 + d, d);

        if (complexity < complexity_min)
        {
            complexity_min = complexity;
            params.K = i + 1;
        }
    }

    return params;
}


}
