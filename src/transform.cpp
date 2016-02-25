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

#include "fgt.hpp"

namespace fgt {

Transform::Transform(const double* source, size_t rows, size_t cols,
                     double bandwidth)
    : m_source(source),
      m_rows_source(rows),
      m_cols(cols),
      m_bandwidth(bandwidth) {}

std::vector<double> Transform::compute(const double* target, size_t rows) {
    std::vector<double> weights(this->rows_source(), 1.0);
    return compute(target, rows, weights.data());
}

std::vector<double> Transform::compute(const double* target, size_t rows,
                                       const double* weights) {
    return compute_impl(target, rows, weights);
}
}
