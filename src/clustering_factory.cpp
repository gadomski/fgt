#include "clustering_factory.hpp"

#include "clustering.hpp"
#include "clustering/gonzalez.hpp"


namespace ifgt
{


ClusteringFactory::ClusteringFactory() {}


Clustering ClusteringFactory::compute(const arma::mat& source, arma::uword k,
                                      double bandwidth, double epsilon)
{
    Clustering clustering = gonzalez_clustering(source, k, bandwidth, epsilon,
                                                false, arma::uword());
    clustering.initialize();
    return clustering;
}
}
