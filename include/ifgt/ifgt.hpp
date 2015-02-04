#include <armadillo>

#include <ifgt/clustering.hpp>


namespace ifgt
{


struct Parameters;


arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon);
arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::vec& q);
arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::vec& q, const Parameters& params);
arma::vec ifgt(const ClusteringUnqPtr& clustering, const arma::mat& Y,
        double h, const arma::vec& q, const Parameters& params);


}
