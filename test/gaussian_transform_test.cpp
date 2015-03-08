#include <fgt/fgt.hpp>

#include <fgt/exceptions.hpp>

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt
{


class MockGaussianTransform : public GaussianTransform {
public:
    using GaussianTransform::GaussianTransform;

private:
    virtual arma::vec compute_impl(const arma::mat& target,
            const arma::vec& weights) const override {
        return arma::zeros<arma::vec>(target.n_rows);
    }

};


TEST(GaussianTransform, IncorrectDimensions) {
    arma::mat source(1, 2);
    arma::mat target(1, 3);
    double bandwidth = 1;
    MockGaussianTransform transform(source, bandwidth);
    EXPECT_THROW(transform.compute(target), dimension_mismatch);
}
}
