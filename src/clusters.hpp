#pragma once

#include <armadillo>

#include "KMlocal.h"


namespace ifgt
{


class Clusters
{
public:

    Clusters(const arma::mat& X, int K);

private:

    KMdata m_data;
    KMfilterCenters m_centers;
    KMterm m_term;
    // TODO do we want to allow the user to pick the algorithm? I don't see
    // why not.
    KMlocalHybrid m_algorithm;
    std::vector<KMctrIdx> m_indices;
    double m_rx;

};


}
