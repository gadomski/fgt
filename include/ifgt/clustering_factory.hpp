#pragma once

#include <armadillo>

#include <ifgt/clustering.hpp>


namespace ifgt
{


class ClusteringFactory
{
public:
    
    ClusteringFactory();

    ClusteringUnqPtr create(const arma::mat& source, arma::uword k, double bandwidth, double epsilon);

};


}
