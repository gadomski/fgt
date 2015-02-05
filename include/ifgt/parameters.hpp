#pragma once

#include <armadillo>


namespace ifgt
{


struct Parameters
{
    int K;
    double r;
};


Parameters choose_parameters(arma::uword d, double bandwidth, double epsilon);
Parameters choose_parameters(arma::uword d, double bandwidth, double epsilon, arma::uword k_limit);


}
