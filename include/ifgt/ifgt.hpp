#include <armadillo>


namespace ifgt
{

class Clustering;
struct Parameters;


arma::vec ifgt(const arma::mat& source, const arma::mat& target,
               double bandwidth, double epsilon);
arma::vec ifgt(const arma::mat& source, const arma::mat& target,
               double bandwidth, double epsilon, const arma::vec& q);
arma::vec ifgt(const arma::mat& source, const arma::mat& target,
               double bandwidth, double epsilon, const arma::vec& q,
               const Parameters& params);
arma::vec ifgt(const Clustering& clustering, const arma::mat& target,
               double bandwidth, const arma::vec& q, const Parameters& params);
}
