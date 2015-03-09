#include <armadillo>

#include <fgt/typedefs.hpp>


namespace fgt {


class GaussianTransform {
public:
    GaussianTransform(const arma::mat& source, double bandwidth);
    virtual ~GaussianTransform();

    arma::vec compute(const arma::mat& target) const;
    arma::vec compute(const arma::mat& target, const arma::vec& weights) const;
    double get_bandwidth() const { return m_bandwidth; }
    arma::uword get_dimensions() const { return m_source.n_cols; }
    const arma::mat& get_source() const { return m_source; }
    arma::uword get_source_n_rows() const { return m_source.n_rows; }

private:
    virtual arma::vec compute_impl(const arma::mat& target,
                                   const arma::vec& weights) const = 0;

    const arma::mat& m_source;
    double m_bandwidth;
};


class Direct : public GaussianTransform {
public:
    using GaussianTransform::GaussianTransform;

private:
    virtual arma::vec compute_impl(const arma::mat& target,
                                   const arma::vec& weights) const override;
};


class DirectTree : public GaussianTransform {
public:
    static const size_t MaxLeafSize = 10;

    DirectTree(const arma::mat& source, double bandwidth, double epsilon);

private:
    virtual arma::vec compute_impl(const arma::mat& target,
                                   const arma::vec& weights) const override;
    // I could set Dimensions to -1 and then determine at runtime.
    template <arma::uword Dimensions>
    arma::vec compute_impl_with_dimensions(const arma::mat& target,
                                           const arma::vec& weights) const;

    double m_epsilon;
    size_t m_max_leaf;
};


class Ifgt : public GaussianTransform {
public:
    struct Parameters {
        arma::uword num_clusters;
        double radius;
    };

    static const arma::uword DefaultNumClustersLimit = 50;
    static const arma::uword MaxNumClusters = 200;
    static const arma::uword NumClusterLimitFactor = 20;
    static const bool DefaultUseDataAdaptive = false;

    Ifgt(const arma::mat& source, double bandwidth, double epsilon,
         int k_limit = DefaultNumClustersLimit);

    static Parameters choose_parameters(arma::uword dimensions,
                                        double bandwidth, double epsilon);
    static Parameters choose_parameters(arma::uword dimensions,
                                        double bandwidth, double epsilon,
                                        arma::uword k_limit);

    optional_arma_uword_t get_clustering_starting_index() const;
    Ifgt& set_clustering_starting_index(arma::uword index);
    Ifgt& use_data_adaptive_truncation(bool data_adaptive);

private:
    virtual arma::vec compute_impl(const arma::mat& target,
                                   const arma::vec& weights) const override;

    double m_epsilon;
    int m_k_limit;
    optional_arma_uword_t m_clustering_starting_index;
    bool m_data_adaptive;
};
}
