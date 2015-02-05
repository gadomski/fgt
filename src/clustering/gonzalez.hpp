#include <ifgt/clustering.hpp>


namespace ifgt
{


Clustering gonzalez_clustering(const arma::mat& source, int K, double bandwidth,
        double epsilon, bool use_starting_idx, arma::uword starting_idx); 


}
