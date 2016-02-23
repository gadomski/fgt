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

namespace fgt {

Clustering cluster(const double* points, size_t rows, size_t cols,
                   size_t nclusters, double epsilon,
                   const std::vector<double>& starting_clusters) {
    std::vector<double> clusters(starting_clusters);
    std::vector<double> temp_clusters(nclusters * cols);
    std::vector<size_t> counts(nclusters);
    std::vector<size_t> labels(rows);
    double error = 0.0;
    double old_error = 0.0;

    do {
        old_error = error;
        error = 0.0;
        temp_clusters.assign(nclusters * cols, 0.0);
        counts.assign(nclusters, 0.0);

        for (size_t i = 0; i < rows; ++i) {
            double min_distance = std::numeric_limits<double>::max();
            for (size_t j = 0; j < nclusters; ++j) {
                double distance = 0;
                for (size_t k = 0; k < cols; ++k) {
                    distance += std::pow(
                        points[i * cols + k] - clusters[j * cols + k], 2);
                }
                if (distance < min_distance) {
                    labels[i] = j;
                    min_distance = distance;
                }
            }

            for (size_t k = 0; k < cols; ++k) {
                temp_clusters[labels[i] * cols + k] += points[i * cols + k];
            }
            ++counts[labels[i]];
            error += min_distance;
        }

        for (size_t j = 0; j < nclusters; ++j) {
            for (size_t k = 0; k < cols; ++k) {
                clusters[j * cols + k] =
                    counts[j] ? temp_clusters[j * cols + k] / counts[j]
                              : temp_clusters[j * cols + k];
            }
        }
    } while (std::abs(error - old_error) > epsilon);

    double max_radius = std::numeric_limits<double>::min();
    std::vector<double> radii(nclusters, std::numeric_limits<double>::min());
    for (size_t i = 0; i < rows; ++i) {
        double distance = 0.0;
        for (size_t k = 0; k < cols; ++k) {
            distance += std::pow(
                points[i * cols + k] - clusters[labels[i] * cols + k], 2);
        }
        distance = std::sqrt(distance);
        if (distance > radii[labels[i]]) {
            radii[labels[i]] = distance;
        }
        if (distance > max_radius) {
            max_radius = distance;
        }
    }

    return {max_radius, labels, clusters, counts, radii};
}
}
