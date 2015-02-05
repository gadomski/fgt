#include <ifgt/clustering_factory.hpp>

#include "clustering/gonzalez.hpp"


namespace ifgt
{


ClusteringFactory::ClusteringFactory()
{}


ClusteringUnqPtr ClusteringFactory::create(const arma::mat& source, arma::uword k,
        double bandwidth, double epsilon)
{
    ClusteringUnqPtr clustering(new Gonzalez(source, k, bandwidth, epsilon, false, arma::uword()));
    clustering->compute();
    return clustering;
}


}
