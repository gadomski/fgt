#pragma once

#include <armadillo>


namespace ifgt
{


class Clusters
{};


Clusters cluster(const arma::mat& X, int K);


}
