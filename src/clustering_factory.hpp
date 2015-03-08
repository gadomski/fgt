#pragma once

#include <armadillo>


namespace ifgt {


class Clustering;


class ClusteringFactory {
public:
    ClusteringFactory();

    Clustering compute(const arma::mat& source, arma::uword k, double bandwidth,
                       double epsilon);
};
}
