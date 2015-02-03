#pragma once

#include <armadillo>


namespace ifgt
{


class Clusters
{
public:

    Clusters(const arma::mat& X, int K);

    inline arma::uvec get_indices() const
    {
        return m_indices;
    }
    inline arma::mat get_centers() const
    {
        return m_centers;
    }
    inline arma::uvec get_num_points() const
    {
        return m_num_points;
    }
    inline arma::vec get_radii() const
    {
        return m_radii;
    }
    inline double get_rx() const
    {
        return m_rx;
    }

    void compute();

private:

    const arma::mat& m_X;
    arma::uvec m_indices;
    arma::mat m_centers;
    arma::uvec m_num_points;
    arma::vec m_radii;
    double m_rx;

};


Clusters compute_clusters(const arma::mat& X, int K);


}
