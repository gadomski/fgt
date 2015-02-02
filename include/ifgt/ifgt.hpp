#include <armadillo>


namespace ifgt
{


struct Parameters;
struct Clusters;


arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon);
arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::mat& q);
arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::mat& q, const Parameters& params);
arma::vec ifgt(const arma::mat& X, const arma::mat& Y, double h, double epsilon,
               const arma::mat& q, const Parameters& params, const Clusters& clusters);


}
