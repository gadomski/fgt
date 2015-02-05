#include "clustering/gonzalez.hpp"


namespace ifgt
{


namespace
{


double ddist(const arma::rowvec& x, const arma::rowvec& y)
{
    return arma::accu(arma::pow(x - y, 2));
}


}



Gonzalez::Gonzalez(const arma::mat& source, int K, double bandwidth, double epsilon,
        bool use_starting_idx, arma::uword starting_idx)
    : Clustering(source, K, bandwidth, epsilon)
    , m_use_starting_idx(use_starting_idx)
    , m_starting_idx(starting_idx)
{}


void Gonzalez::cluster()
{
    arma::uword N = get_indices().n_rows;
    arma::uword K = get_centers().n_rows;
    arma::uvec centers(K);
    arma::uvec cprev(N);
    arma::uvec cnext(N);
    arma::uvec far2c(K);
    arma::vec dist(N);

    arma::uword nc;
    if (m_starting_idx)
    {
        nc = m_starting_idx;
    }
    else
    {
        std::random_device rd;
        nc = rd() % N;
    }
    centers(0) = nc;

    for (arma::uword i = 0; i < N; ++i)
    {
        dist(i) = (i == nc) ? 0.0 : ddist(get_source_row(i), get_source_row(nc));
        cnext(i) = i + 1;
        cprev(i) = i - 1;
    }

    cnext(N - 1) = 0;
    cprev(0) = N - 1;

    nc = std::max_element(dist.begin(), dist.end()) - dist.begin();
    far2c(0) = nc;
    set_radius(0, dist(nc));

    for (int i = 1; i < K; ++i)
    {
        nc = far2c(get_radius_idxmax(i));

        centers(i) = nc;
        set_radius(i, 0.0);
        dist(nc) = 0.0;
        set_index(nc, i);
        far2c(i) = nc;

        cnext(cprev(nc)) = cnext(nc);
        cprev(cnext(nc)) = cprev(nc);
        cnext(nc) = nc;
        cprev(nc) = nc;

        for (int j = 0; j < i; ++j)
        {
            arma::uword ct_j = centers(j);
            double dc2cq = ddist(get_source_row(ct_j), get_source_row(nc)) / 4;
            if (dc2cq < get_radius(j))
            {
                set_radius(j, 0.0);
                far2c(j) = ct_j;
                arma::uword k = cnext(ct_j);
                while (k != ct_j)
                {
                    arma::uword nextk = cnext(k);
                    double dist2c_k = dist(k);
                    if (dc2cq < dist2c_k)
                    {
                        double dd = ddist(get_source_row(k), get_source_row(nc));
                        if (dd < dist2c_k)
                        {
                            dist(k) = dd;
                            set_index(k, i);
                            if (get_radius(i) < dd)
                            {
                                set_radius(i, dd);
                                far2c(i) = k;
                            }
                            cnext(cprev(k)) = nextk;
                            cprev(nextk) = cprev(k);
                            cnext(k) = cnext(nc);
                            cprev(cnext(nc)) = k;
                            cnext(nc) = k;
                            cprev(k) = nc;
                        }
                        else if (get_radius(j) < dist2c_k)
                        {
                            set_radius(j, dist2c_k);
                            far2c(j) = k;
                        }
                    }
                    else if (get_radius(j) < dist2c_k)
                    {
                        set_radius(j, dist2c_k);
                        far2c(j) = k;
                    }
                    k = nextk;
                }
            }
        }
    }

    set_radii(arma::sqrt(get_radii()));
    set_rx(get_max_radius());

    arma::mat center_coordinates(get_centers());
    for (arma::uword i = 0; i < N; ++i)
    {
        increment_num_points(get_index(i));
        center_coordinates.row(get_index(i)) += get_source_row(i);
    }
    set_centers(center_coordinates / arma::repmat(get_num_points(), 1, get_d()));
}


}
