#include <ifgt/ifgt.hpp>

#include "choose_parameters.hpp"
#include "clusters.hpp"


namespace ifgt
{


arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon)
{
    arma::vec q = arma::zeros<arma::vec>(X.n_rows);
    return ifgt(X, Y, h, epsilon, q);
}


arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::mat& q)
{
    Parameters params = choose_parameters(X.n_cols, h, epsilon);
    return ifgt(X, Y, h, epsilon, q, params);
}


arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::mat& q, const Parameters& params)
{
    Clusters clusters = cluster(X, params.k);
    return ifgt(X, Y, h, epsilon, q, params, clusters);
}


arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::mat& q, const Parameters& params, const Clusters& clusters)
{
    // TODO check X.n_cols == Y.n_cols
    throw std::runtime_error("Not implemented");
    return arma::vec();
}


}
