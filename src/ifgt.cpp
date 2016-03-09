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
#include <vector>

#include "cluster.hpp"
#include "fgt.hpp"
#include "ifgt.hpp"

namespace fgt {
namespace {

// TODO make this configurable
const size_t TRUNCATION_NUMBER_UL = 200;

Matrix::Index nchoosek(Matrix::Index n, Matrix::Index k) {
    Matrix::Index n_k = n - k;
    if (k < n_k) {
        k = n_k;
        n_k = n - k;
    }
    Matrix::Index nchsk = 1;
    for (Matrix::Index i = 1; i <= n_k; ++i) {
        nchsk *= ++k;
        nchsk /= i;
    }
    return nchsk;
}
}

Vector ifgt(const MatrixRef source, const MatrixRef target, double bandwidth,
            double epsilon) {
    return Ifgt(source, bandwidth, epsilon).compute(target);
}

Vector ifgt(const MatrixRef source, const MatrixRef target, double bandwidth,
            double epsilon, const VectorRef weights) {
    return Ifgt(source, bandwidth, epsilon).compute(target, weights);
}

IfgtParameters ifgt_choose_parameters(Matrix::Index cols, double bandwidth,
                                      double epsilon,
                                      Matrix::Index max_num_clusters,
                                      Matrix::Index truncation_number_ul) {
    double h2 = bandwidth * bandwidth;
    double radius = bandwidth * std::sqrt(std::log(1.0 / epsilon));
    double complexity_min = std::numeric_limits<double>::max();
    Matrix::Index nclusters = 0;

    for (Matrix::Index i = 0; i < max_num_clusters; ++i) {
        double rx = std::pow(double(i + 1), -1.0 / double(cols));
        double rx2 = rx * rx;
        double n = std::min(double(i + 1), std::pow(radius / rx, double(cols)));
        double error = std::numeric_limits<double>::max();
        double temp = 1.0;
        Matrix::Index p = 0;
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

Matrix::Index
ifgt_choose_truncation_number(Matrix::Index cols, double bandwidth,
                              double epsilon, double rx,
                              Matrix::Index truncation_number_ul) {
    double h2 = bandwidth * bandwidth;
    double rx2 = rx * rx;
    double r = std::min(std::sqrt(cols),
                        bandwidth * std::sqrt(std::log(1.0 / epsilon)));
    double error = std::numeric_limits<double>::max();
    Matrix::Index p = 0;
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

Ifgt::Ifgt(const MatrixRef source, double bandwidth, double epsilon)
    : Transform(source, bandwidth),
      m_epsilon(epsilon),
      m_nclusters(0),
      m_clustering(),
      m_truncation_number(0),
      m_p_max_total(0),
      m_constant_series() {
    // TODO max num clusters should be configurable
    Matrix::Index max_num_clusters(
        Matrix::Index(std::round(0.2 * 100 / bandwidth)));
    IfgtParameters params =
        ifgt_choose_parameters(source.cols(), bandwidth, epsilon,
                               max_num_clusters, TRUNCATION_NUMBER_UL);
    if (params.nclusters == 0) {
        throw ifgt_no_clusters();
    }
    m_nclusters = params.nclusters;
    // TODO make the clustering constructor do the work?
    m_clustering.reset(new Clustering(cluster(source, m_nclusters, epsilon)));
    m_truncation_number = ifgt_choose_truncation_number(
        source.cols(), bandwidth, epsilon, m_clustering->max_radius,
        TRUNCATION_NUMBER_UL);
    m_p_max_total =
        nchoosek(m_truncation_number - 1 + source.cols(), source.cols());
    m_constant_series = compute_constant_series();
    m_ry_square.resize(m_nclusters);
    for (Matrix::Index j = 0; j < m_nclusters; ++j) {
        double ry = params.cutoff_radius + m_clustering->radii[j];
        m_ry_square[j] = ry * ry;
    }
}

Ifgt::~Ifgt() {}

Vector Ifgt::compute_monomials(const VectorRef d) const {
    auto cols = this->source().cols();
    std::vector<Matrix::Index> heads(unsigned(cols), 0);
    Vector monomials = Vector::Ones(p_max_total());
    for (Matrix::Index k = 1, t = 1, tail = 1; k < m_truncation_number;
         ++k, tail = t) {
        for (Matrix::Index i = 0; i < cols; ++i) {
            Matrix::Index head = heads[unsigned(i)];
            heads[unsigned(i)] = t;
            for (Matrix::Index j = head; j < tail; ++j, ++t) {
                monomials[t] = d[i] * monomials[j];
            }
        }
    }
    return monomials;
}

Vector Ifgt::compute_constant_series() const {
    auto cols = this->source().cols();
    std::vector<Matrix::Index> heads(unsigned(cols + 1), 0);
    heads[unsigned(cols)] = std::numeric_limits<Matrix::Index>::max();
    std::vector<Matrix::Index> cinds(unsigned(p_max_total()), 0);
    Vector monomials = Vector::Ones(p_max_total());

    for (Matrix::Index k = 1, t = 1, tail = 1; k < m_truncation_number;
         ++k, tail = t) {
        for (Matrix::Index i = 0; i < cols; ++i) {
            Matrix::Index head = heads[unsigned(i)];
            heads[unsigned(i)] = t;
            for (Matrix::Index j = head; j < tail; ++j, ++t) {
                cinds[unsigned(t)] =
                    (j < heads[unsigned(i) + 1]) ? cinds[unsigned(j)] + 1 : 1;
                monomials[t] = 2.0 * monomials[j];
                monomials[t] /= double(cinds[unsigned(t)]);
            }
        }
    }
    return monomials;
}

Vector Ifgt::compute_impl(const MatrixRef target,
                          const VectorRef weights) const {
    auto source = this->source();
    auto rows_source = source.rows();
    auto rows_target = target.rows();
    auto cols = source.cols();
    auto bandwidth = this->bandwidth();
    auto nclusters = this->nclusters();
    auto p_max_total = this->p_max_total();

    double h2 = bandwidth * bandwidth;

    Matrix C = Matrix::Zero(nclusters, p_max_total);
    for (Matrix::Index i = 0; i < rows_source; ++i) {
        double distance = 0.0;
        Vector dx = Vector::Zero(cols);
        for (Matrix::Index k = 0; k < cols; ++k) {
            double delta = source(i, k) -
                           m_clustering->clusters(m_clustering->indices[i], k);
            distance += delta * delta;
            dx[k] = delta / bandwidth;
        }

        auto monomials = compute_monomials(dx);
        double f = weights[i] * std::exp(-distance / h2);
        for (Matrix::Index alpha = 0; alpha < p_max_total; ++alpha) {
            C(m_clustering->indices[i], alpha) += f * monomials[alpha];
        }
    }

#pragma omp parallel for
    for (Matrix::Index j = 0; j < nclusters; ++j) {
        for (Matrix::Index alpha = 0; alpha < p_max_total; ++alpha) {
            C(j, alpha) *= m_constant_series[alpha];
        }
    }

    Vector G = Vector::Zero(rows_target);
#pragma omp parallel for
    for (Matrix::Index i = 0; i < rows_target; ++i) {
        for (Matrix::Index j = 0; j < nclusters; ++j) {
            double distance = 0.0;
            Vector dy = Vector::Zero(cols);
            for (Matrix::Index k = 0; k < cols; ++k) {
                double delta = target(i, k) - m_clustering->clusters(j, k);
                distance += delta * delta;
                if (distance > m_ry_square[j]) {
                    break;
                }
                dy[k] = delta / bandwidth;
            }
            if (distance <= m_ry_square[j]) {
                auto monomials = compute_monomials(dy);
                double g = std::exp(-distance / h2);
                G[i] += C.row(j) * g * monomials;
            }
        }
    }

    return G;
}
}
