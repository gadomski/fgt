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
///
/// This header includes both a functional interface and a class-based interface
/// to the fgt library.

#pragma once

#include <cstddef>
#include <memory>
#include <vector>

/// Top-level namespace for all things fgt.
namespace fgt {

/// Returns the version of the fgt library as a string.
const char* version();

/// Returns true if the library was compiled with OpenMP support.
bool with_openmp();

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

/// Abstract base class for all supported variants of the Gauss transform.
///
/// Some flavors of transform can pre-compute some data, e.g. the `DirectTree`
/// can pre-compute the KD-tree.
/// This pre-computation allows reuse of those data structure for multiple runs
/// of the transform, potentially with different target data sets.
class Transform {
public:
    /// Constructs a new transform that can be re-used with different targets.
    Transform(const double* source, size_t rows, size_t cols, double bandwidth);

    /// Returns the pointer to the source dataset.
    const double* source() const { return m_source; }
    /// Returns the number of rows in the source dataset.
    size_t rows_source() const { return m_rows_source; }
    /// Returns the number of columns.
    size_t cols() const { return m_cols; }
    /// Returns the bandwidth of the transform.
    double bandwidth() const { return m_bandwidth; }

    /// Computes the Gauss transform for the given target dataset.
    std::vector<double> compute(const double* target, size_t rows);
    /// Computes the Gauss transform with the given weights.
    std::vector<double> compute(const double* target, size_t rows,
                                const double* weights);

private:
    virtual std::vector<double> compute_impl(const double* target, size_t rows,
                                             const double* weights) const = 0;

    const double* m_source;
    size_t m_rows_source;
    size_t m_cols;
    double m_bandwidth;
};

/// Direct Gauss transform.
class Direct : public Transform {
public:
    /// Creates a new direct transform.
    Direct(const double* source, size_t rows, size_t cols, double bandwidth);

private:
    virtual std::vector<double> compute_impl(const double* target, size_t rows,
                                             const double* weights) const;
};

/// Direct Gauss transform using a KD-tree truncation.
class DirectTree : public Transform {
public:
    /// Creates a new direct tree transform.
    ///
    /// This constructor pre-computes the KD-tree, so subsequent calls to
    /// `compute()` will re-use the same tree.
    DirectTree(const double* source, size_t rows, size_t cols, double bandwidth,
               double epsilon);

    /// Destroys a DirectTree.
    ///
    /// Required because of the unique pointer to a incomplete class.
    ~DirectTree();

    /// Returns the error tolerance value.
    double epsilon() const { return m_epsilon; }

private:
    struct NanoflannTree;

    virtual std::vector<double> compute_impl(const double* target, size_t rows,
                                             const double* weights) const;

    double m_epsilon;
    std::unique_ptr<NanoflannTree> m_tree;
};

/// Opaque k-means clustering structure.
struct Clustering;

/// Improved Fast Gauss Transform.
class Ifgt : public Transform {
public:
    /// Creates a new Ifgt.
    ///
    /// This constructor will precompute some values, including the clusters and
    /// monomials, in hopes of speeding up subsequent runs.
    Ifgt(const double* source, size_t rows, size_t cols, double bandwidth,
         double epsilon);

    /// Destroys this transform.
    ///
    /// Required because PIMPL.
    ~Ifgt();

    /// Returns the error tolerance value.
    double epsilon() const { return m_epsilon; }
    /// Returns the number of clusters.
    size_t nclusters() const { return m_nclusters; }
    /// Returns the truncation number.
    size_t truncation_number() const { return m_truncation_number; }
    /// Returns the length of each monomial.
    size_t p_max_total() const { return m_p_max_total; }

private:
    virtual std::vector<double> compute_impl(const double* target, size_t rows,
                                             const double* weights) const;
    std::vector<double> compute_monomials(const std::vector<double>& d) const;
    std::vector<double> compute_constant_series() const;

    double m_epsilon;
    size_t m_nclusters;
    std::unique_ptr<Clustering> m_clustering;
    size_t m_truncation_number;
    size_t m_p_max_total;
    std::vector<double> m_constant_series;
    std::vector<double> m_ry_square;
};
}
