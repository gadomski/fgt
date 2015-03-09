#include <fgt/fgt.hpp>

#include <fgt/exceptions.hpp>


namespace fgt {


GaussianTransform::GaussianTransform(const arma::mat& source, double bandwidth)
    : m_source(source), m_bandwidth(bandwidth) {}


GaussianTransform::~GaussianTransform() {}


arma::vec GaussianTransform::compute(const arma::mat& target) const {
    return compute(target, arma::ones<arma::vec>(get_source_n_rows()));
}

arma::vec GaussianTransform::compute(const arma::mat& target,
                                     const arma::vec& weights) const {
    if (m_source.n_cols != target.n_cols) {
        std::stringstream ss;
        ss << "Dimentionality of source and target do not match ("
           << m_source.n_cols << " vs " << target.n_cols << ")";
        throw dimension_mismatch(ss.str());
    }
    if (m_source.n_rows != weights.n_rows) {
        std::stringstream ss;
        ss << "Source and weights do not have the same number of rows ("
           << m_source.n_rows << " vs " << weights.n_rows << ")";
        throw dimension_mismatch(ss.str());
    }
    return compute_impl(target, weights);
}
}
