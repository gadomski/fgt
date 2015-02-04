#pragma once

#include <armadillo>

#include <ifgt/parameters.hpp>


namespace ifgt
{


Parameters choose_parameters(arma::uword d, double h, double epsilon);
Parameters choose_parameters(arma::uword d, double h, double epsilon, arma::uword k_limit);


}
