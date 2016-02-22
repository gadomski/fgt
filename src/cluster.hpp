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

#include <cstddef>
#include <vector>

namespace fgt {

/// The results from k-means clustering.
struct Clustering {
    /// The maximum cluster radius.
    double max_radius;
    /// The cluster membership ids for each points.
    std::vector<size_t> indices;
    /// The centers of each cluster.
    std::vector<double> clusters;
    /// The number of points in each cluster.
    std::vector<size_t> npoints;
    /// The radius of each cluster.
    std::vector<double> radii;
};

/// Runs k-means clustering on a set of points.
Clustering cluster(const double* points, size_t rows, size_t cols,
                   size_t nclusters, double epsilon);

/// Runs k-means clustering, specifying the starting cluster centers.
Clustering cluster(const double* points, size_t rows, size_t cols,
                   size_t nclusters, double epsilon,
                   const std::vector<double>& starting_clusters);
}
