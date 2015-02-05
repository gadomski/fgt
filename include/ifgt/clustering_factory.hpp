#pragma once

#include <armadillo>

#include <ifgt/clustering.hpp>


namespace ifgt
{


class ClusteringFactory
{
public:
    
    ClusteringFactory();

    Clustering compute(const arma::mat& source, arma::uword k, double bandwidth, double epsilon);

};


}
