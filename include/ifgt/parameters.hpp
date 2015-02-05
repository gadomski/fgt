#pragma once

#include <armadillo>


namespace ifgt
{


struct Parameters
{
    int K;
    double r;
};


Parameters choose_parameters(arma::uword d, double h, double epsilon);
Parameters choose_parameters(arma::uword d, double h, double epsilon, arma::uword k_limit);


}
