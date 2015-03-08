#include <fgt/fgt.hpp>

#include "clustering/gonzalez.hpp"
#include "clustering.hpp"
#include "config.hpp"
#include "parameters.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt {


TEST(Ifgt, ReferenceImplementation) {
    arma::mat source, target;
    source.load(test_data_path("X.csv"));
    target.load(test_data_path("Y.csv"));
    double bandwidth = 0.4;
    double epsilon = 1e-3;
    arma::vec weights = arma::ones<arma::vec>(source.n_rows);
    int k_limit = 50;

    Ifgt ifgt(source, bandwidth, epsilon, k_limit);
    ifgt.set_clustering_starting_index(2);
    arma::vec g = ifgt.compute(target, weights);

    EXPECT_EQ(g.n_rows, 5000);
    EXPECT_DOUBLE_EQ(346.36423735983732, arma::min(g));
    EXPECT_DOUBLE_EQ(2190.6088281193192, arma::max(g));
}
}
