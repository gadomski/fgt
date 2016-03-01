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

#include <Eigen/Core>

/// Top-level namespace for all things fgt.
namespace fgt {

/// Convenience typedef for our type of matrix.
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
    Matrix;

/// Convenience typedef for a reference to our type of Eigen::Matrix.
///
/// We accept references to our functions in case the argument isn't exactly a
/// Eigen::Matrix, for whatever reason.
typedef Eigen::Ref<const Matrix> MatrixRef;

/// Convenience typedef for our flavor of Eigen::Vector.
typedef Eigen::VectorXd Vector;

/// Convenience typedef for a reference to a Vector.
typedef Eigen::Ref<const Vector> VectorRef;

/// Returns the version of the fgt library as a string.
const char* version();

/// Returns true if the library was compiled with OpenMP support.
bool with_openmp();

/// Computes the direct Gauss transform with equal weights.
Vector direct(const MatrixRef source, const MatrixRef target, double bandwidth);

/// Computes the direct Gauss transform with provided weights.
Vector direct(const MatrixRef source, const MatrixRef target, double bandwidth,
              const VectorRef weights);

/// Computes the direct Gauss transform using a kd-tree.
Vector direct_tree(const MatrixRef source, const MatrixRef target,
                   double bandwidth, double epsilon);

/// Computes the direct Gauss transform using a kd-tree and weights.
Vector direct_tree(const MatrixRef source, const MatrixRef target,
                   double bandwidth, double epsilon, const VectorRef weights);

/// Computes the Improved Fast Gauss Transform.
Vector ifgt(const MatrixRef source, const MatrixRef target, double bandwidth,
            double epsilon);

/// Computes the Improved Fast Gauss Transform with the provided weights.
Vector ifgt(const MatrixRef source, const MatrixRef target, double bandwidth,
            double epsilon, const VectorRef weights);

/// Abstract base class for all supported variants of the Gauss transform.
///
/// Some flavors of transform can pre-compute some data, e.g. the `DirectTree`
/// can pre-compute the KD-tree.
/// This pre-computation allows reuse of those data structure for multiple runs
/// of the transform, potentially with different target data sets.
class Transform {
public:
    /// Constructs a new transform that can be re-used with different targets.
    Transform(const MatrixRef source, double bandwidth);

    /// Returns the pointer to the source dataset.
    const MatrixRef source() const { return m_source; }
    /// Returns the bandwidth of the transform.
    double bandwidth() const { return m_bandwidth; }

    /// Computes the Gauss transform for the given target dataset.
    Vector compute(const MatrixRef target);
    /// Computes the Gauss transform with the given weights.
    Vector compute(const MatrixRef target, const VectorRef weights);

private:
    virtual Vector compute_impl(const MatrixRef target,
                                const VectorRef weights) const = 0;

    const MatrixRef m_source;
    double m_bandwidth;
};

/// Direct Gauss transform.
class Direct : public Transform {
public:
    /// Creates a new direct transform.
    Direct(const MatrixRef source, double bandwidth);

private:
    virtual Vector compute_impl(const MatrixRef target,
                                const VectorRef weights) const;
};

/// Direct Gauss transform using a KD-tree truncation.
class DirectTree : public Transform {
public:
    /// Creates a new direct tree transform.
    ///
    /// This constructor pre-computes the KD-tree, so subsequent calls to
    /// `compute()` will re-use the same tree.
    DirectTree(const MatrixRef source, double bandwidth, double epsilon);

    /// Destroys a DirectTree.
    ///
    /// Required because of the unique pointer to a incomplete class.
    ~DirectTree();

    /// Returns the error tolerance value.
    double epsilon() const { return m_epsilon; }

private:
    struct NanoflannTree;

    virtual Vector compute_impl(const MatrixRef target,
                                const VectorRef weights) const;

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
    Ifgt(const MatrixRef source, double bandwidth, double epsilon);

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
    virtual Vector compute_impl(const MatrixRef target,
                                const VectorRef weights) const;
    Vector compute_monomials(const VectorRef d) const;
    Vector compute_constant_series() const;

    double m_epsilon;
    size_t m_nclusters;
    std::unique_ptr<Clustering> m_clustering;
    size_t m_truncation_number;
    size_t m_p_max_total;
    Vector m_constant_series;
    Vector m_ry_square;
};
}
