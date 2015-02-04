#include <ifgt/ifgt.hpp>

#include <ifgt/choose_parameters.hpp>
#include <ifgt/clustering_factory.hpp>

#include "monomials.hpp"


namespace ifgt
{


arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon)
{
    arma::vec q = arma::ones<arma::vec>(X.n_rows);
    return ifgt(X, Y, h, epsilon, q);
}


arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::vec& q)
{
    Parameters params = choose_parameters(X.n_cols, h, epsilon);
    return ifgt(X, Y, h, epsilon, q, params);
}


arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::vec& q, const Parameters& params)
{
    ClusteringFactory factory;
    ClusteringUnqPtr clustering = factory.create(X, params.K, h, epsilon);
    return ifgt(clustering, Y, h, q, params);
}


arma::vec ifgt(const ClusteringUnqPtr& clustering, const arma::mat& Y, double h,
               const arma::vec& q, const Parameters& params)
{
    // TODO check X.n_cols == Y.n_cols
    arma::vec G(Y.n_rows);
    arma::vec ry2 = arma::pow(params.r + clustering->get_radii(), 2);
    double h2 = h * h;
    arma::mat C = clustering->find_C(q);
    for (arma::uword j = 0; j < Y.n_rows; ++j)
    {
        G(j) = 0.0;
        for (arma::uword k = 0; k < params.K; ++k)
        {
            arma::rowvec dy = Y.row(j) - clustering->get_centers().row(k);
            double distance2 = arma::accu(arma::pow(dy, 2));
            if (distance2 <= ry2(k))
            {
                double g = std::exp(-distance2 / h2);
                G(j) += arma::accu(C.row(k) % compute_monomials(dy / h, clustering->get_p_max()) * g);
            }
        }
    }
    return G;
}


}
