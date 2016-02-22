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

#include <cassert>
#include <cmath>

#include "fgt.hpp"

namespace fgt {

std::vector<double> direct(const double* source, size_t rows_source,
                           const double* target, size_t rows_target,
                           size_t cols, double bandwidth) {
    std::vector<double> weights(rows_source, 1.0);
    return direct(source, rows_source, target, rows_target, cols, bandwidth,
                  weights.data());
}

std::vector<double> direct(const double* source, size_t rows_source,
                           const double* target, size_t rows_target,
                           size_t cols, double bandwidth,
                           const double* weights) {
    double h2 = bandwidth * bandwidth;
    std::vector<double> g(rows_target, 0.0);
#pragma omp parallel for
    for (size_t j = 0; j < rows_target; ++j) {
        for (size_t i = 0; i < rows_source; ++i) {
            double distance = 0.0;
            for (size_t k = 0; k < cols; ++k) {
                double temp = source[cols * i + k] - target[cols * j + k];
                distance += temp * temp;
            }
            g[j] += weights[i] * std::exp(-distance / h2);
        }
    }
    return g;
}
}
