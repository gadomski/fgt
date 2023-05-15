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

Transform::Transform(const MatrixRef source, double bandwidth)
    : m_source(source), m_bandwidth(bandwidth) {}

Vector Transform::compute(const MatrixRef target) {
    Vector weights = Vector::Ones(this->source().rows());
    return compute(target, weights);
}

Vector Transform::compute(const MatrixRef target, const VectorRef weights) {
    return compute_impl(target, weights);
}

Matrix Transform::matrix_compute(const MatrixRef target) {
  return matrix_compute_impl(target);
}
}
