#include <fgt/fgt.hpp>


namespace fgt {


GaussianTransform::GaussianTransform(const arma::mat& source, double bandwidth)
    : m_source(source), m_bandwidth(bandwidth) {}


GaussianTransform::~GaussianTransform() {}


arma::vec GaussianTransform::compute(const arma::mat& target) const {
    return this->compute(target, arma::ones<arma::vec>(get_source_n_rows()));
}
}
