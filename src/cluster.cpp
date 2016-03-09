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
#include <random>

#include "cluster.hpp"
#include "fgt.hpp"

namespace fgt {

Matrix pick_cluster_centers(const MatrixRef points, Matrix::Index nclusters) {
    std::default_random_engine generator;
    std::uniform_int_distribution<Matrix::Index> distribution(0, points.rows() -
                                                                     1);
    Matrix::Index cols = points.cols();
    Matrix clusters(nclusters, cols);
    for (Matrix::Index j = 0; j < nclusters; ++j) {
        Matrix::Index index = distribution(generator);
        for (Matrix::Index k = 0; k < cols; ++k) {
            clusters(j, k) = points(index, k);
        }
    }
    return clusters;
}

Clustering cluster(const MatrixRef points, Matrix::Index nclusters,
                   double epsilon) {
    Matrix clusters = pick_cluster_centers(points, nclusters);
    return cluster(points, nclusters, epsilon, clusters);
}
}
