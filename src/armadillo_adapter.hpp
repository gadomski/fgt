#pragma once

#include <armadillo>
#include <nanoflann.hpp>


namespace fgt {


class ArmadilloAdaptor {
public:
    ArmadilloAdaptor(const arma::mat& data) : m_data(data) {}

    const arma::mat& derived() const { return m_data; };
    arma::uword kdtree_get_point_count() const { return m_data.n_rows; }
    double kdtree_distance(const double* p1, const arma::uword idx_p2,
                           size_t size) const {
        return arma::accu(
            arma::pow(m_data.row(idx_p2) - arma::rowvec(p1, size), 2));
    }
    double kdtree_get_pt(const arma::uword idx, int dim) const {
        return m_data(idx, dim);
    }
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& bounding_box) const {
        return false;
    }

private:
    const arma::mat& m_data;
};


typedef nanoflann::L2_Simple_Adaptor<double, ArmadilloAdaptor>
    L2_Simple_ArmadilloAdaptor;
}
