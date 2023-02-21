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

Direct::Direct(const MatrixRef source, double bandwidth)
    : Transform(source, bandwidth) {}

Vector Direct::compute_impl(const MatrixRef target,
                            const VectorRef weights) const {
    double h2 = bandwidth() * bandwidth();
    MatrixRef source = this->source();
    auto rows_source = source.rows();
    auto rows_target = target.rows();
    Vector g = Vector::Zero(rows_target);
#pragma omp parallel for
    for (Matrix::Index j = 0; j < rows_target; ++j) {
        for (Matrix::Index i = 0; i < rows_source; ++i) {
            double distance =
                (source.row(i) - target.row(j)).array().pow(2).sum();
            g[j] += weights[i] * std::exp(-distance / h2);
        }
    }
    return g;
}

Matrix Direct::matrix_compute_impl(const MatrixRef target) const {
    double h2 = bandwidth() * bandwidth();
    MatrixRef source = this->source();
    auto rows_source = source.rows();
    auto rows_target = target.rows();
    Matrix g = Matrix::Zero(rows_target, rows_source);
#pragma omp parallel for
    for (Matrix::Index j = 0; j < rows_target; ++j) {
        for (Matrix::Index i = 0; i < rows_source; ++i) {
            double distance =
                (source.row(i) - target.row(j)).array().pow(2).sum();
            g(j, i) = std::exp(-distance / h2);
        }
    }
    return g;
}

Vector direct(const MatrixRef source, const MatrixRef target,
              double bandwidth) {
    return Direct(source, bandwidth).compute(target);
}

Vector direct(const MatrixRef source, const MatrixRef target, double bandwidth,
              const VectorRef weights) {
    return Direct(source, bandwidth).compute(target, weights);
}

Matrix mat_direct(const MatrixRef source, const MatrixRef target, double bandwidth) {
  return Direct(source, bandwidth).matrix_compute(target);
}
}
