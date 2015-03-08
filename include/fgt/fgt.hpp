#include <armadillo>

#include <fgt/typedefs.hpp>


namespace fgt {


class GaussianTransform {
public:
    GaussianTransform(const arma::mat& source, double bandwidth);
    virtual ~GaussianTransform();

    arma::vec compute(const arma::mat& target) const;
    virtual arma::vec compute(const arma::mat& target,
                              const arma::vec& weights) const = 0;
    double get_bandwidth() const { return m_bandwidth; }
    const arma::mat& get_source() const { return m_source; }
    arma::uword get_source_n_rows() const { return m_source.n_rows; }

private:
    const arma::mat& m_source;
    double m_bandwidth;
};


class Direct : public GaussianTransform {
public:
    Direct(const arma::mat& source, double bandwidth);

    using GaussianTransform::compute;
    virtual arma::vec compute(const arma::mat& target,
                              const arma::vec& weights) const override;
};


class Ifgt : public GaussianTransform {
public:
    static const int DefaultKLimit = 50;

    Ifgt(const arma::mat& source, double bandwidth, double epsilon,
         int k_limit = DefaultKLimit);

    using GaussianTransform::compute;
    virtual arma::vec compute(const arma::mat& target,
                              const arma::vec& weights) const override;
    optional_arma_uword_t get_clustering_starting_index() const;
    Ifgt& set_clustering_starting_index(arma::uword index);

private:
    double m_epsilon;
    int m_k_limit;
    optional_arma_uword_t m_clustering_starting_index;
};
}
