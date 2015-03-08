#include <fgt/fgt.hpp>

#include "clustering.hpp"
#include "monomials.hpp"
#include "nchoosek.hpp"

#include <algorithm>
#include <cmath>
#include <limits>


namespace fgt {

Ifgt::Parameters Ifgt::choose_parameters(arma::uword d, double bandwidth, double epsilon) {
    return choose_parameters(d, bandwidth, epsilon,
                             std::round(NumClusterLimitFactor / bandwidth));
}


Ifgt::Parameters Ifgt::choose_parameters(arma::uword d, double bandwidth, double epsilon,
                             arma::uword k_limit) {
    Parameters params;
    double R = std::sqrt(d);
    double h2 = bandwidth * bandwidth;
    double complexity_min = std::numeric_limits<double>::max();
    double rx = 0.0;

    params.radius = std::min(R, bandwidth * std::sqrt(std::log(1 / epsilon)));
    params.num_clusters = 1;

    for (arma::uword i = 0; i < k_limit; ++i) {
        rx = std::pow(double(i + 1), -1.0 / double(d));
        double rx2 = rx * rx;
        double n = std::min(double(i + 1), std::pow(params.radius / rx, double(d)));
        double error = 1;
        double temp = 1;
        arma::uword p = 0;

        while ((error > epsilon) and (p <= MaxNumClusters)) {
            ++p;
            double b =
                std::min((rx + std::sqrt(rx2 + 2 * p * h2)) / 2, rx + params.radius);
            double c = rx - b;
            temp *= 2 * rx * b / h2 / p;
            error = temp * std::exp(-c * c / h2);
        }
        double complexity =
            i + 1 + std::log(double(i + 1)) + (n + 1) * nchoosek(p - 1 + d, d);

        if (complexity < complexity_min) {
            complexity_min = complexity;
            params.num_clusters = i + 1;
        }
    }

    return params;
}


Ifgt::Ifgt(const arma::mat& source, double bandwidth, double epsilon,
           int k_limit)
    : GaussianTransform(source, bandwidth),
      m_epsilon(epsilon),
      m_k_limit(k_limit),
      m_clustering_starting_index(std::make_pair(false, 0)) {}


optional_arma_uword_t Ifgt::get_clustering_starting_index() const {
    return m_clustering_starting_index;
}


Ifgt& Ifgt::set_clustering_starting_index(arma::uword index) {
    m_clustering_starting_index = std::make_pair(true, index);
    return *this;
}


arma::vec Ifgt::compute(const arma::mat& target,
                        const arma::vec& weights) const {
    const arma::mat& source = get_source();
    double bandwidth = get_bandwidth();
    Parameters params = choose_parameters(source.n_cols, bandwidth, m_epsilon);
    GonzalezClustering clustering(source, params.num_clusters, bandwidth, m_epsilon,
                                  get_clustering_starting_index());
    // TODO check source.n_cols == target.n_cols
    arma::vec G(target.n_rows);
    arma::vec ry2 = arma::pow(params.radius + clustering.get_radii(), 2);
    double h2 = bandwidth * bandwidth;
    arma::mat C = clustering.compute_C(weights);
    for (arma::uword j = 0; j < target.n_rows; ++j) {
        G(j) = 0.0;
        for (arma::uword k = 0; k < params.num_clusters; ++k) {
            arma::rowvec dy = target.row(j) - clustering.get_centers().row(k);
            double distance2 = arma::accu(arma::pow(dy, 2));
            if (distance2 <= ry2(k)) {
                double g = std::exp(-distance2 / h2);
                G(j) += arma::accu(
                    C.row(k) %
                    compute_monomials(dy / bandwidth, clustering.get_p_max()) *
                    g);
            }
        }
    }
    return G;
}
}
