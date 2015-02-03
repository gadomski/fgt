#pragma once

#include <memory>

#include <armadillo>


namespace ifgt
{


class Clustering
{
public:

    Clustering(const arma::mat& X, const arma::vec& q, int K, double h, double epsilon);

    virtual ~Clustering();

    inline arma::rowvec get_X_row(arma::uword i) const
    {
        return m_X.row(i);
    }
    inline const arma::uvec& get_indices() const
    {
        return m_indices;
    }
    inline arma::uword get_index(arma::uword i) const
    {
        return m_indices(i);
    }
    inline const arma::mat& get_centers() const
    {
        return m_centers;
    }
    inline const arma::uvec& get_num_points() const
    {
        return m_num_points;
    }
    inline const arma::vec& get_radii() const
    {
        return m_radii;
    }
    inline double get_radius(arma::uword i) const
    {
        return m_radii(i);
    }
    inline arma::uword get_radius_idxmax(arma::uword i) const
    {
        return std::max_element(m_radii.begin(), m_radii.begin() + i) - m_radii.begin();
    }
    inline double get_max_radius() const
    {
        return *std::max_element(m_radii.begin(), m_radii.end());
    }
    inline double get_rx() const
    {
        return m_rx;
    }
    inline double get_p_max() const
    {
        return m_p_max;
    }
    inline const arma::mat& get_C() const
    {
        return m_C;
    }
    inline arma::uword get_d() const
    {
        return m_X.n_cols;
    }

    inline void set_radius(arma::uword i, double r)
    {
        m_radii(i) = r;
    }
    inline void set_radii(arma::vec radii)
    {
        m_radii = radii;
    }
    inline void set_index(arma::uword i, arma::uword index)
    {
        m_indices(i) = index;
    }
    inline void set_rx(double rx)
    {
        m_rx = rx;
    }
    inline void increment_num_points(arma::uword i)
    {
        ++m_num_points(i);
    }
    inline void set_centers(const arma::mat& centers)
    {
        m_centers = centers;
    }

    void compute();

private:

    virtual void cluster() = 0;
    void compute_C();

    const arma::mat& m_X;
    const arma::vec& m_q;
    arma::uvec m_indices;
    arma::mat m_centers;
    arma::uvec m_num_points;
    arma::vec m_radii;
    double m_rx;
    double m_h;
    double m_epsilon;
    arma::uword m_p_max;
    arma::mat m_C;

};


typedef std::unique_ptr<Clustering> ClusteringUnqPtr;


}
