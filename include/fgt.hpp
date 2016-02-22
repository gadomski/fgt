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

/// \file fgt.hpp
/// \brief The header file for the fgt library.

#pragma once

#include <cstddef>
#include <vector>

namespace fgt {

/// Returns the version of the fgt library as a string.
const char* version();

/// Returns the result of `git describe` for this library's version.
const char* git_describe();

/// Computes the direct Gauss transform with equal weights.
std::vector<double> direct(const double* source, size_t rows_source,
                           const double* target, size_t rows_target,
                           size_t cols, double bandwidth);

/// Computes the direct Gauss transform with provided weights.
std::vector<double> direct(const double* source, size_t rows_source,
                           const double* target, size_t rows_target,
                           size_t cols, double bandwidth,
                           const double* weights);

/// Computes the direct Gauss transform using a kd-tree.
std::vector<double> direct_tree(const double* source, size_t rows_source,
                                const double* target, size_t rows_target,
                                size_t cols, double bandwidth, double epsilon);

/// Computes the direct Gauss transform using a kd-tree and weights.
std::vector<double> direct_tree(const double* source, size_t rows_source,
                                const double* target, size_t rows_target,
                                size_t cols, double bandwidth, double epsilon,
                                const double* weights);

/// Computes the Improved Fast Gauss Transform.
std::vector<double> ifgt(const double* source, size_t rows_source,
                         const double* target, size_t rows_target, size_t cols,
                         double bandwidth, double epsilon);

/// Computes the Improved Fast Gauss Transform with the provided weights.
std::vector<double> ifgt(const double* source, size_t rows_source,
                         const double* target, size_t rows_target, size_t cols,
                         double bandwidth, double epsilon,
                         const double* weights);
}
