#include <fgt/fgt.hpp>

#include "armadillo_adapter.hpp"

#include <armadillo>
#include <nanoflann.hpp>

#include <cmath>
#include <sstream>
#include <utility>
#include <vector>


namespace fgt {


DirectTree::DirectTree(const arma::mat& source, double bandwidth,
                       double epsilon)
    : GaussianTransform(source, bandwidth),
      m_epsilon(epsilon),
      m_max_leaf(MaxLeafSize) {}


arma::vec DirectTree::compute_impl(const arma::mat& target,
                                   const arma::vec& weights) const {
    double bandwidth2 = get_bandwidth() * get_bandwidth();
    double cutoff_radius = get_bandwidth() * std::sqrt(std::log(1 / m_epsilon));
    double cutoff_radius2 = cutoff_radius * cutoff_radius;
    arma::vec g = arma::zeros<arma::vec>(target.n_rows);
    typedef nanoflann::KDTreeSingleIndexAdaptor<
        L2_Simple_ArmadilloAdaptor, ArmadilloAdaptor, -1, arma::uword> tree_t;
    ArmadilloAdaptor adapter(get_source());
    tree_t tree(get_dimensions(), adapter,
                nanoflann::KDTreeSingleIndexAdaptorParams(m_max_leaf));
    tree.buildIndex();

    std::vector<double> point(get_dimensions());
    std::vector<std::pair<arma::uword, double>> indices_distances;
    indices_distances.reserve(get_source_n_rows());
    nanoflann::SearchParams search_params;
    search_params.sorted = false;
    for (int j = 0; j < target.n_rows; ++j) {
        point =
            std::move(arma::conv_to<std::vector<double>>::from(target.row(j)));
        size_t num_points_found = tree.radiusSearch(
            point.data(), cutoff_radius2, indices_distances, search_params);
        for (size_t i = 0; i < num_points_found; ++i) {
            auto entry = indices_distances[i];
            g(j) += weights(entry.first) * std::exp(-entry.second / bandwidth2);
        }
    }
    return g;
}
}
