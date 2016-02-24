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

Direct::Direct(const double* source, size_t rows, size_t cols, double bandwidth)
    : Transform(source, rows, cols, bandwidth) {}

std::vector<double> Direct::compute_impl(const double* target,
                                         size_t rows_target,
                                         const double* weights) const {
    double h2 = bandwidth() * bandwidth();
    const double* source = this->source();
    size_t rows_source = this->rows_source();
    size_t cols = this->cols();
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

std::vector<double> direct(const double* source, size_t rows_source,
                           const double* target, size_t rows_target,
                           size_t cols, double bandwidth) {
    return Direct(source, rows_source, cols, bandwidth)
        .compute(target, rows_target);
}

std::vector<double> direct(const double* source, size_t rows_source,
                           const double* target, size_t rows_target,
                           size_t cols, double bandwidth,
                           const double* weights) {
    return Direct(source, rows_source, cols, bandwidth)
        .compute(target, rows_target, weights);
}
}
