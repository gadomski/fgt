#include <fgt/fgt.hpp>

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt {


class Trial : public ::testing::TestWithParam<double> {
protected:
    static const arma::uword NumSourceRows = 1000;
    static const arma::uword NumTargetRows = 1000;
    static const arma::uword NumCols = 3;
    constexpr static const double Epsilon = 1e-2;

    Trial()
        : ::testing::TestWithParam<double>(),
          m_source(arma::randu<arma::mat>(NumSourceRows, NumCols)),
          m_target(arma::randu<arma::mat>(NumTargetRows, NumCols)),
          m_weights(arma::randu<arma::vec>(NumSourceRows)),
          m_weights_sum(arma::accu(m_weights)),
          m_epsilon(Epsilon) {}

    arma::mat m_source, m_target;
    arma::vec m_weights;
    double m_weights_sum;
    double m_epsilon;
};


TEST_P(Trial, AllMethods) {
    Direct direct(m_source, GetParam());
    arma::vec g_direct = direct.compute(m_target, m_weights);

    Ifgt ifgt(m_source, GetParam(), m_epsilon);
    arma::vec g_ifgt = ifgt.compute(m_source, m_weights);

    double error_ifgt = arma::max(arma::abs(g_ifgt - g_direct) / m_weights_sum);
    EXPECT_GT(m_epsilon, error_ifgt);
}


INSTANTIATE_TEST_CASE_P(Trial, Trial,
                        ::testing::Values(0.01, 0.02, 0.04, 0.08, 0.16, 0.32,
                                          0.64));
}
