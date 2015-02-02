#pragma once

#include <armadillo>

#include "KMlocal.h"


namespace ifgt
{


struct Clusters
{
    KMfilterCenters centers;
};


Clusters cluster(const arma::mat& X, int K);


}
