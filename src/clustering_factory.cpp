#include <ifgt/clustering_factory.hpp>

#include "clustering/gonzalez.hpp"


namespace ifgt
{


ClusteringFactory::ClusteringFactory()
{}


ClusteringUnqPtr ClusteringFactory::create(const arma::mat& X, arma::uword k,
        double h, double epsilon)
{
    ClusteringUnqPtr clustering(new Gonzalez(X, k, h, epsilon, false, arma::uword()));
    clustering->compute();
    return clustering;
}


}
