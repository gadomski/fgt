// fgt, C++ library for Fast Gauss Transforms
// Copyright (C) 2015 Peter J. Gadomski <pete.gadomski@gmail.com>
//
// This library is free software; you can redistribute it and/or modify it
// under the terms of the GNU Lesser General Public License as published by the
// Free Software Foundation; either version 2.1 of the License, or (at your
// option) any later version.
//
// This library is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation,
// Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301 USA

/// \file fgt.hpp
/// \brief Header file for the fgt C++ library.
///
/// Include this file to declare everything you need to run fgt in C++.

#pragma once

#include <chrono>
#include <memory>
#include <stdexcept>

#include <armadillo>

/// Top-level namespace for all things fgt.
namespace fgt {

/// Used when you *might* have a armadillo size value.
typedef std::pair<bool, arma::uword> optional_arma_uword_t;

/// Top-level error class for all fgt errors.
class fgt_error : public std::runtime_error {
public:
    using std::runtime_error::runtime_error;
};

/// Two matrices do not have compatible dimensions.
class dimension_mismatch : public fgt_error {
public:
    using fgt_error::fgt_error;
};

/// A generic gauss transform.
///
/// This is an abstract superclass for all flavors of gauss transform.
class GaussTransform {
public:
    /// The result of a Gauss transform.
    ///
    /// Along with the transform values themselves, this result contains
    /// information about the actual run.
    struct Result {
        arma::vec data;
        std::chrono::microseconds runtime;
    };

    /// Creates a new Gauss transform using the source points and the given
    /// bandwidth.
    GaussTransform(const arma::mat& source, double bandwidth);

    /// Destroys a Gauss transform and frees all associated memory.
    virtual ~GaussTransform();

    /// Computes the Gauss transform between the source and target matrices.
    Result compute(const arma::mat& target) const;

    /// Computes the Gauss transform between the source and target matrices,
    /// applying the given weights to the target matrices.
    Result compute(const arma::mat& target, const arma::vec& weights) const;

    /// Returns this transform's bandwidth.
    double get_bandwidth() const { return m_bandwidth; }

    /// Returns the number of columns in the source matrix.
    arma::uword get_dimensions() const { return m_source.n_cols; }

    /// Returns a constant reference to the source matrix.
    const arma::mat& get_source() const { return m_source; }

    /// Returns the number of rows in the source matrix.
    arma::uword get_source_n_rows() const { return m_source.n_rows; }

private:
    virtual arma::vec compute_impl(const arma::mat& target,
                                   const arma::vec& weights) const = 0;

    const arma::mat& m_source;
    double m_bandwidth;
};

/// The direct Gauss transform.
///
/// This is the naive implementation of the Gauss transform, which is O^2 for
/// the number of points in the source and target matrices.
/// This method is, however, exact, and can therefore be used as a reference.
class Direct : public GaussTransform {
public:
    using GaussTransform::GaussTransform;

private:
    virtual arma::vec compute_impl(const arma::mat& target,
                                   const arma::vec& weights) const override;
};

/// The direct Gauss transform, but simplified using clustering.
///
/// Since far-away points do not contribute much to the total gauss transform at
/// a target point, this method uses clustering to only compute the transform
/// for "close" points.
class DirectTree : public GaussTransform {
public:
    /// The maximum number of leaves in the clustering.
    static const size_t MaxLeafSize = 10;

    /// Creates a new transform with the given error tolerance, epsilon.
    DirectTree(const arma::mat& source, double bandwidth, double epsilon);

private:
    virtual arma::vec compute_impl(const arma::mat& target,
                                   const arma::vec& weights) const override;

    double m_epsilon;
    size_t m_max_leaf;
};

/// The Improved fast Gauss Transform.
///
/// The improved transform is an epsilon-exact approximation of the Gauss
/// transform.
/// In general, it performs better when the source and target matrices are of
/// higher dimensions.
class Ifgt : public GaussTransform {
public:
    /// Reusable parameter object.
    struct Parameters {
        /// The number of clusters that should be used for the IFGT.
        arma::uword num_clusters;
        /// The maximum radius that we should see when we do the clustering.
        double radius;
    };

    /// The maximum number of clusters that we can use.
    static const arma::uword MaxNumClusters = 200;
    /// A factor used to choose the k_limit for choosing parameters, if none is
    /// provided.
    static const arma::uword NumClusterLimitFactor = 20;
    /// Default value for data-adaptive IFGT.
    static const bool DefaultUseDataAdaptive = false;

    /// Creates a new IFGT object, with the given error tolerance (epsilon).
    Ifgt(const arma::mat& source, double bandwidth, double epsilon);
    /// Creates a new IFGT object and specifies whether it should be data
    /// adaptive.
    Ifgt(const arma::mat& source, double bandwidth, double epsilon,
         bool data_adaptive);

    /// Choose the correct IFGT parameters for the given dimensionality,
    /// bandwidth, and error tolerance.
    static Parameters choose_parameters(arma::uword dimensions,
                                        double bandwidth, double epsilon);
    /// Choose the correct IFGT parameters for the given dimensionality,
    /// bandwidth, error tolerance, and maximum number of clusters.
    static Parameters choose_parameters(arma::uword dimensions,
                                        double bandwidth, double epsilon,
                                        arma::uword k_limit);

    /// Returns an optional clustering starting index.
    ///
    /// If this is not set, the clustering will start at a random spot, making
    /// things hard to test exactly.
    optional_arma_uword_t get_clustering_starting_index() const;
    /// Sets the clustering starting index.
    Ifgt& set_clustering_starting_index(arma::uword index);
    /// Should we use data adaptive truncation?
    ///
    /// Data adaptive truncation requires a bit of extra front-end work, but can
    /// reduce computational complexity during the algorithm run itself.
    Ifgt& use_data_adaptive_truncation(bool data_adaptive);

private:
    virtual arma::vec compute_impl(const arma::mat& target,
                                   const arma::vec& weights) const override;

    double m_epsilon;
    optional_arma_uword_t m_clustering_starting_index;
    bool m_data_adaptive;
};

/// Unique pointer to a `GaussTransform`.
typedef std::unique_ptr<GaussTransform> GaussTransformUnqPtr;

/// Choose The Best™ Gauss transform for these data.
///
/// The best transform is a bit arbitrary, and at this point we just have some
/// hard-coded cutoff points.
/// You are encouraged to explicitly pick the transform of your choice if you
/// have a better sense of what you want to use.
GaussTransformUnqPtr choose_gauss_transform(const arma::mat& source,
                                            double bandwidth, double epsilon);
}
