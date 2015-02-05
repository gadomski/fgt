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



Clustering gonzalez_clustering(const arma::mat& source, int K, double bandwidth, double epsilon,
        bool use_starting_idx, arma::uword starting_idx)
{
    Clustering clustering(source, K, bandwidth, epsilon);

    arma::uword N = source.n_rows;
    arma::uvec centers(K);
    arma::uvec cprev(N);
    arma::uvec cnext(N);
    arma::uvec far2c(K);
    arma::vec dist(N);

    arma::uword nc;
    if (use_starting_idx)
    {
        nc = starting_idx;
    }
    else
    {
        std::random_device rd;
        nc = rd() % N;
    }
    centers(0) = nc;

    for (arma::uword i = 0; i < N; ++i)
    {
        dist(i) = (i == nc) ? 0.0 : ddist(source.row(i), source.row(nc));
        cnext(i) = i + 1;
        cprev(i) = i - 1;
    }

    cnext(N - 1) = 0;
    cprev(0) = N - 1;

    nc = std::max_element(dist.begin(), dist.end()) - dist.begin();
    far2c(0) = nc;
    clustering.set_radius(0, dist(nc));

    for (int i = 1; i < K; ++i)
    {
        nc = far2c(clustering.get_radius_idxmax(i));

        centers(i) = nc;
        clustering.set_radius(i, 0.0);
        dist(nc) = 0.0;
        clustering.set_index(nc, i);
        far2c(i) = nc;

        cnext(cprev(nc)) = cnext(nc);
        cprev(cnext(nc)) = cprev(nc);
        cnext(nc) = nc;
        cprev(nc) = nc;

        for (int j = 0; j < i; ++j)
        {
            arma::uword ct_j = centers(j);
            double dc2cq = ddist(source.row(ct_j), source.row(nc)) / 4;
            if (dc2cq < clustering.get_radius(j))
            {
                clustering.set_radius(j, 0.0);
                far2c(j) = ct_j;
                arma::uword k = cnext(ct_j);
                while (k != ct_j)
                {
                    arma::uword nextk = cnext(k);
                    double dist2c_k = dist(k);
                    if (dc2cq < dist2c_k)
                    {
                        double dd = ddist(source.row(k), source.row(nc));
                        if (dd < dist2c_k)
                        {
                            dist(k) = dd;
                            clustering.set_index(k, i);
                            if (clustering.get_radius(i) < dd)
                            {
                                clustering.set_radius(i, dd);
                                far2c(i) = k;
                            }
                            cnext(cprev(k)) = nextk;
                            cprev(nextk) = cprev(k);
                            cnext(k) = cnext(nc);
                            cprev(cnext(nc)) = k;
                            cnext(nc) = k;
                            cprev(k) = nc;
                        }
                        else if (clustering.get_radius(j) < dist2c_k)
                        {
                            clustering.set_radius(j, dist2c_k);
                            far2c(j) = k;
                        }
                    }
                    else if (clustering.get_radius(j) < dist2c_k)
                    {
                        clustering.set_radius(j, dist2c_k);
                        far2c(j) = k;
                    }
                    k = nextk;
                }
            }
        }
    }

    clustering.set_radii(arma::sqrt(clustering.get_radii()));
    clustering.set_rx(clustering.get_max_radius());

    arma::mat center_coordinates(clustering.get_centers());
    for (arma::uword i = 0; i < N; ++i)
    {
        clustering.increment_num_points(clustering.get_index(i));
        center_coordinates.row(clustering.get_index(i)) += source.row(i);
    }
    clustering.set_centers(center_coordinates /
            arma::repmat(clustering.get_num_points(), 1, clustering.get_d()));

    return clustering;
}


}
