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

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable : 4820) // paddding
#pragma warning(disable : 4626) // assignment implicitly deleted
#pragma warning(disable : 5027) // move assignment implicitly deleted
#pragma warning(disable : 4127) // conditional expression is constant
#pragma warning(disable : 4365) // signed/unsigned
#endif
#include "nanoflann.hpp"
#ifdef _MSC_VER
#pragma warning(pop)
#endif

#include "fgt.hpp"

namespace fgt {
namespace {

struct MatrixAdaptor {
    size_t kdtree_get_point_count() const { return size_t(m_rows); }
    double kdtree_distance(const double* p1, const Matrix::Index idx_p2,
                           Matrix::Index) const {
        double distance = 0.0;
        for (Matrix::Index k = 0; k < m_cols; ++k) {
            double temp = p1[k] - m_data[idx_p2 * m_cols + k];
            distance += temp * temp;
        }
        return distance;
    }
    double kdtree_get_pt(const Matrix::Index idx, Matrix::Index dim) const {
        return m_data[m_cols * idx + dim];
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX&) const {
        return false;
    }

    const double* m_data;
    Matrix::Index m_rows;
    Matrix::Index m_cols;
};
}

Vector direct_tree(const MatrixRef source, const MatrixRef target,
                   double bandwidth, double epsilon) {
    return DirectTree(source, bandwidth, epsilon).compute(target);
}

Vector direct_tree(const MatrixRef source, const MatrixRef target,
                   double bandwidth, double epsilon, const VectorRef weights) {
    return DirectTree(source, bandwidth, epsilon).compute(target, weights);
}

struct DirectTree::NanoflannTree {
    typedef nanoflann::KDTreeSingleIndexAdaptor<
        nanoflann::L2_Simple_Adaptor<double, MatrixAdaptor>, MatrixAdaptor>
        tree_t;

    NanoflannTree(const MatrixRef source)
        : matrix_adaptor({source.data(), source.rows(), source.cols()}),
          tree(int(source.cols()), matrix_adaptor) {}

    NanoflannTree(const NanoflannTree&) = delete;
    NanoflannTree& operator=(const NanoflannTree&) = delete;
    NanoflannTree& operator=(NanoflannTree&&) = delete;

    MatrixAdaptor matrix_adaptor;
    tree_t tree;
};

DirectTree::DirectTree(const MatrixRef source, double bandwidth, double epsilon)
    : Transform(source, bandwidth),
      m_epsilon(epsilon),
      m_tree(new DirectTree::NanoflannTree(source)) {
    m_tree->tree.buildIndex();
}

DirectTree::~DirectTree() {}

Vector DirectTree::compute_impl(const MatrixRef target,
                                const VectorRef weights) const {
    double h2 = bandwidth() * bandwidth();
    double cutoff_radius = bandwidth() * std::sqrt(std::log(1.0 / epsilon()));
    double r2 = cutoff_radius * cutoff_radius;
    Matrix::Index rows_source = this->source().rows();
    Matrix::Index rows_target = target.rows();
    Vector g = Vector::Zero(rows_target);
    Matrix::Index cols = this->source().cols();

    nanoflann::SearchParams params;
    params.sorted = false;

#pragma omp parallel for
    for (Matrix::Index j = 0; j < rows_target; ++j) {
        std::vector<std::pair<size_t, double>> indices_distances;
        indices_distances.reserve(unsigned(rows_source));
        size_t nfound = m_tree->tree.radiusSearch(&target.data()[j * cols], r2,
                                                  indices_distances, params);
        for (size_t i = 0; i < nfound; ++i) {
            auto entry = indices_distances[i];
            g[j] += weights[signed(entry.first)] * std::exp(-entry.second / h2);
        }
    }
    return g;
}
}
