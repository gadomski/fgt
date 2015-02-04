#pragma once

#include <armadillo>

#include <ifgt/clustering.hpp>


namespace ifgt
{


class ClusteringFactory
{
public:
    
    ClusteringFactory();

    ClusteringUnqPtr create(const arma::mat& X, arma::uword k, double h, double epsilon);

};


}
