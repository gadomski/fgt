#include <gtest/gtest.h>
#include <ifgt/ifgt.hpp>

#include <armadillo>

#include "choose_parameters.hpp"
#include "clustering/gonzalez.hpp"

#include "config.hpp"


namespace ifgt
{


TEST(Ifgt, ReferenceImplementation)
{
    arma::mat X;
    X.load(test_data_path("X.csv"));
    arma::mat Y;
    Y.load(test_data_path("Y.csv"));

    double h = 0.4;
    double epsilon = 1e-3;
    arma::vec q = arma::ones<arma::vec>(X.n_rows);

    Parameters params = choose_parameters(X.n_cols, h, epsilon, 50);
    ClusteringUnqPtr clustering(new Gonzalez(X, q, params.K, h, epsilon, true, 2));
    clustering->compute();

    arma::vec G = ifgt(clustering, Y, h, params);

    EXPECT_EQ(G.n_rows, 5000);
    EXPECT_NEAR(359.4559, arma::min(G), 0.0001);
    EXPECT_NEAR(2.2214e3, arma::max(G), 0.1);
}


}
