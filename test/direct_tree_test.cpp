#include <fgt/fgt.hpp>

#include "config.hpp"

#include <armadillo>
#include <gtest/gtest.h>


namespace fgt
{


TEST(DirectTree, ReferenceImplementation)
{
    arma::mat source, target;
    source.load(test_data_path("X.csv"));
    target.load(test_data_path("Y.csv"));
    double bandwidth = 0.4;
    double epsilon = 1e-3;
    arma::vec weights = arma::ones<arma::vec>(source.n_rows);

    DirectTree direct_tree(source, bandwidth, epsilon);
    arma::vec g = direct_tree.compute(target, weights);

    EXPECT_EQ(5000, g.n_rows);
    EXPECT_DOUBLE_EQ(2071.8710052017841, g(0));
}


}
