// fgt — fast Gauss transforms
// Copyright (C) 2016 Peter J. Gadomski <pete.gadomski@gmail.com>
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA

#include <cmath>
#include <limits>

#include "cluster.hpp"
#include "fgt.hpp"
#include "ifgt.hpp"

namespace fgt {
namespace {

const size_t TRUNCATION_NUMBER_UL = 200;

int nchoosek(int n, int k) {
    int n_k = n - k;
    if (k < n_k) {
        k = n_k;
        n_k = n - k;
    }
    int nchsk = 1;
    for (int i = 1; i <= n_k; ++i) {
        nchsk *= ++k;
        nchsk /= i;
    }
    return nchsk;
}
}

std::vector<double> ifgt(const double* source, size_t rows_source,
                         const double* target, size_t rows_target, size_t cols,
                         double bandwidth, double epsilon) {
    return Ifgt(source, rows_source, cols, bandwidth, epsilon)
        .compute(target, rows_target);
}

std::vector<double> ifgt(const double* source, size_t rows_source,
                         const double* target, size_t rows_target, size_t cols,
                         double bandwidth, double epsilon,
                         const double* weights) {
    return Ifgt(source, rows_source, cols, bandwidth, epsilon)
        .compute(target, rows_target, weights);
}

IfgtParameters ifgt_choose_parameters(size_t cols, double bandwidth,
                                      double epsilon, size_t max_num_clusters,
                                      size_t truncation_number_ul) {
    double h2 = bandwidth * bandwidth;
    double radius = std::min(std::sqrt(cols),
                             bandwidth * std::sqrt(std::log(1.0 / epsilon)));
    double complexity_min = std::numeric_limits<double>::max();
    size_t nclusters = 0;

    for (size_t i = 0; i < max_num_clusters; ++i) {
        double rx = std::pow(double(i + 1), -1.0 / double(cols));
        double rx2 = rx * rx;
        double n = std::min(double(i + 1), std::pow(radius / rx, double(cols)));
        double error = std::numeric_limits<double>::max();
        double temp = 1.0;
        size_t p = 0;
        while ((error > epsilon) && (p <= truncation_number_ul)) {
            ++p;
            double b =
                std::min((rx + std::sqrt(rx2 + 2.0 * double(p) * h2)) / 2.0,
                         rx + radius);
            double c = rx - b;
            temp *= 2 * rx * b / h2 / double(p);
            error = temp * std::exp(-(c * c) / h2);
        }
        double complexity = i + 1 + std::log(double(i + 1)) +
                            (n + 1) * nchoosek(p - 1 + cols, cols);
        if (complexity < complexity_min) {
            complexity_min = complexity;
            nclusters = i + 1;
        }
    }
    return {nclusters, radius};
}

size_t ifgt_choose_truncation_number(size_t cols, double bandwidth,
                                     double epsilon, double rx,
                                     size_t truncation_number_ul) {
    double h2 = bandwidth * bandwidth;
    double rx2 = rx * rx;
    double r = std::min(std::sqrt(cols),
                        bandwidth * std::sqrt(std::log(1.0 / epsilon)));
    double error = std::numeric_limits<double>::max();
    size_t p = 0;
    double temp = 1.0;
    while ((error > epsilon) && (p <= truncation_number_ul)) {
        ++p;
        double b =
            std::min((rx + std::sqrt(rx2 + 2 * double(p) * h2)) / 2.0, rx + r);
        double c = rx - b;
        temp *= 2 * rx * b / h2 / double(p);
        error = temp * std::exp(-(c * c) / h2);
    }
    return p;
}

Ifgt::Ifgt(const double* source, size_t rows, size_t cols, double bandwidth,
           double epsilon)
    : Transform(source, rows, cols, bandwidth),
      m_epsilon(epsilon),
      m_nclusters(0),
      m_clustering(),
      m_truncation_number(0),
      m_p_max_total(0),
      m_constant_series() {
    size_t max_num_clusters(std::round(0.2 * 100 / bandwidth));
    IfgtParameters params = ifgt_choose_parameters(
        cols, bandwidth, epsilon, max_num_clusters, TRUNCATION_NUMBER_UL);
    m_nclusters = params.nclusters;
    m_clustering.reset(
        new Clustering(cluster(source, rows, cols, m_nclusters, epsilon)));
    m_truncation_number = ifgt_choose_truncation_number(
        cols, bandwidth, epsilon, m_clustering->max_radius,
        TRUNCATION_NUMBER_UL);
    m_p_max_total = nchoosek(m_truncation_number - 1 + cols, cols);
    m_constant_series = compute_constant_series();
    m_ry_square.resize(m_nclusters);
    for (size_t j = 0; j < m_nclusters; ++j) {
        double ry = params.cutoff_radius + m_clustering->radii[j];
        m_ry_square[j] = ry * ry;
    }
}

Ifgt::~Ifgt() {}

std::vector<double>
Ifgt::compute_monomials(const std::vector<double>& d) const {
    auto cols = this->cols();
    std::vector<size_t> heads(cols, 0);
    std::vector<double> monomials(p_max_total(), 1.0);
    for (size_t k = 1, t = 1, tail = 1; k < m_truncation_number;
         ++k, tail = t) {
        for (size_t i = 0; i < cols; ++i) {
            size_t head = heads[i];
            heads[i] = t;
            for (size_t j = head; j < tail; ++j, ++t) {
                monomials[t] = d[i] * monomials[j];
            }
        }
    }
    return monomials;
}

std::vector<double> Ifgt::compute_constant_series() const {
    auto cols = this->cols();
    std::vector<size_t> heads(cols + 1, 0);
    heads[cols] = std::numeric_limits<size_t>::max();
    std::vector<size_t> cinds(p_max_total(), 0);
    std::vector<double> monomials(p_max_total(), 1.0);

    for (size_t k = 1, t = 1, tail = 1; k < m_truncation_number;
         ++k, tail = t) {
        for (size_t i = 0; i < cols; ++i) {
            size_t head = heads[i];
            heads[i] = t;
            for (size_t j = head; j < tail; ++j, ++t) {
                cinds[t] = (j < heads[i + 1]) ? cinds[j] + 1 : 1;
                monomials[t] = 2.0 * monomials[j];
                monomials[t] /= double(cinds[t]);
            }
        }
    }
    return monomials;
}

std::vector<double> Ifgt::compute_impl(const double* target, size_t rows_target,
                                       const double* weights) const {
    auto source = this->source();
    auto rows_source = this->rows_source();
    auto cols = this->cols();
    auto bandwidth = this->bandwidth();
    auto nclusters = this->nclusters();
    auto p_max_total = this->p_max_total();

    double h2 = bandwidth * bandwidth;

    std::vector<double> C(nclusters * p_max_total, 0.0);
    for (size_t i = 0; i < rows_source; ++i) {
        double distance = 0.0;
        std::vector<double> dx(cols, 0.0);
        for (size_t k = 0; k < cols; ++k) {
            double delta =
                source[i * cols + k] -
                m_clustering->clusters[m_clustering->indices[i] * cols + k];
            distance += delta * delta;
            dx[k] = delta / bandwidth;
        }

        auto monomials = compute_monomials(dx);
        double f = weights[i] * std::exp(-distance / h2);
        for (size_t alpha = 0; alpha < p_max_total; ++alpha) {
            C[m_clustering->indices[i] * p_max_total + alpha] +=
                f * monomials[alpha];
        }
    }

#pragma omp parallel for
    for (size_t j = 0; j < nclusters; ++j) {
        for (size_t alpha = 0; alpha < p_max_total; ++alpha) {
            C[j * p_max_total + alpha] *= m_constant_series[alpha];
        }
    }

    std::vector<double> G(rows_target, 0.0);
#pragma omp parallel for
    for (size_t i = 0; i < rows_target; ++i) {
        for (size_t j = 0; j < nclusters; ++j) {
            double distance = 0.0;
            std::vector<double> dy(cols, 0.0);
            for (size_t k = 0; k < cols; ++k) {
                double delta =
                    target[i * cols + k] - m_clustering->clusters[j * cols + k];
                distance += delta * delta;
                if (distance > m_ry_square[j]) {
                    break;
                }
                dy[k] = delta / bandwidth;
            }
            if (distance <= m_ry_square[j]) {
                auto monomials = compute_monomials(dy);
                double g = std::exp(-distance / h2);
                for (size_t alpha = 0; alpha < p_max_total; ++alpha) {
                    G[i] += C[j * p_max_total + alpha] * g * monomials[alpha];
                }
            }
        }
    }

    return G;
}
}
