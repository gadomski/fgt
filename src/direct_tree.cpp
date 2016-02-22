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

#include "nanoflann.hpp"

#include "fgt.hpp"

namespace fgt {
namespace {

struct MatrixAdaptor {
    size_t kdtree_get_point_count() const { return m_rows; }
    double kdtree_distance(const double* p1, const size_t idx_p2,
                           size_t) const {
        double distance = 0.0;
        for (size_t k = 0; k < m_cols; ++k) {
            double temp = p1[k] - m_data[idx_p2 * m_cols + k];
            distance += temp * temp;
        }
        return distance;
    }
    double kdtree_get_pt(const size_t idx, int dim) const {
        return m_data[m_cols * idx + dim];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& bb) const {
        return false;
    }

    const double* m_data;
    size_t m_rows;
    size_t m_cols;
};
}

std::vector<double> direct_tree(const double* source, size_t rows_source,
                                const double* target, size_t rows_target,
                                size_t cols, double bandwidth, double epsilon) {
    std::vector<double> weights(rows_source, 1.0);
    return direct_tree(source, rows_source, target, rows_target, cols,
                       bandwidth, epsilon, weights.data());
}

std::vector<double> direct_tree(const double* source, size_t rows_source,
                                const double* target, size_t rows_target,
                                size_t cols, double bandwidth, double epsilon,
                                const double* weights) {
    double h2 = bandwidth * bandwidth;
    double cutoff_radius = bandwidth * std::sqrt(std::log(1.0 / epsilon));
    double r2 = cutoff_radius * cutoff_radius;
    std::vector<double> g(rows_target, 0.0);

    typedef nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, MatrixAdaptor>, MatrixAdaptor>
        tree_t;
    MatrixAdaptor adapted_source = {source, rows_source, cols};
    tree_t index(cols, adapted_source);
    index.buildIndex();

    nanoflann::SearchParams params;
    params.sorted = false;

#pragma omp parallel for
    for (size_t j = 0; j < rows_target; ++j) {
        std::vector<std::pair<size_t, double>> indices_distances;
        indices_distances.reserve(rows_source);
        size_t nfound = index.radiusSearch(&target[j * cols], r2,
                                           indices_distances, params);
        for (size_t i = 0; i < nfound; ++i) {
            auto entry = indices_distances[i];
            g[j] += weights[i] * std::exp(-entry.second / h2);
        }
    }
    return g;
}
}
