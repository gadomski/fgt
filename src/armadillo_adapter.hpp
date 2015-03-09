#pragma once

#include <armadillo>
#include <nanoflann.hpp>


namespace fgt {


class ArmadilloAdaptor {
public:
    ArmadilloAdaptor(const arma::mat& data) : m_data(data) {}

    arma::uword kdtree_get_point_count() const { return m_data.n_rows; }
    double kdtree_distance(const double* p1, const arma::uword idx_p2,
                           size_t size) const {
        double distance = 0;
        // Direct pointer access used in a (probably misguided) attempt at speed
        const double* point = m_data.memptr() + idx_p2;
        for (size_t i = 0; i < size; ++i) {
            double d = p1[i] - point[i * m_data.n_rows];
            distance += d * d;
        }
        return distance;
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
