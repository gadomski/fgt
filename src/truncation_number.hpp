#pragma once

#include <armadillo>


namespace fgt {


static const arma::uword MaxTruncationNumber = 300;


arma::uword choose_truncation_number(
    double distance2, double cutoff_radius, double epsilon, double bandwidth,
    arma::uword max_truncation_number = MaxTruncationNumber);
}
