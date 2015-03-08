#pragma once

#include <armadillo>


namespace ifgt
{


class Clustering
{
public:
    Clustering(const arma::mat& source, int K, double bandwidth,
               double epsilon);

    arma::rowvec get_source_row(arma::uword i) const { return m_source.row(i); }
    const arma::uvec& get_indices() const { return m_indices; }
    arma::uword get_index(arma::uword i) const { return m_indices(i); }
    const arma::mat& get_centers() const { return m_centers; }
    const arma::uvec& get_num_points() const { return m_num_points; }
    const arma::vec& get_radii() const { return m_radii; }
    double get_radius(arma::uword i) const { return m_radii(i); }
    arma::uword get_radius_idxmax(arma::uword i) const
    {
        return std::max_element(m_radii.begin(), m_radii.begin() + i) -
               m_radii.begin();
    }
    double get_max_radius() const
    {
        return *std::max_element(m_radii.begin(), m_radii.end());
    }
    double get_rx() const { return m_rx; }
    double get_p_max() const { return m_p_max; }
    arma::uword get_d() const { return m_source.n_cols; }

    void set_radius(arma::uword i, double r) { m_radii(i) = r; }
    void set_radii(arma::vec radii) { m_radii = radii; }
    void set_index(arma::uword i, arma::uword index) { m_indices(i) = index; }
    void set_rx(double rx) { m_rx = rx; }
    void increment_num_points(arma::uword i) { ++m_num_points(i); }
    void set_centers(const arma::mat& centers) { m_centers = centers; }

    void initialize();
    arma::mat compute_C(const arma::vec& q) const;

private:
    const arma::mat& m_source;
    arma::uvec m_indices;
    arma::mat m_centers;
    arma::uvec m_num_points;
    arma::vec m_radii;
    double m_rx;
    double m_bandwidth;
    double m_epsilon;
    arma::uword m_p_max;
    arma::rowvec m_constant_series;
    bool m_is_initialized;
};
}
