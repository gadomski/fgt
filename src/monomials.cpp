#include "monomials.hpp"

#include "p_max_total.hpp"


namespace ifgt
{


arma::rowvec compute_monomials(arma::rowvec dx, arma::uword p_max)
{
    arma::vec heads = arma::zeros<arma::vec>(dx.n_cols); 
    arma::rowvec monomials = arma::ones<arma::rowvec>(get_p_max_total(dx.n_cols, p_max));

    for (arma::uword k = 1, t = 1, tail = 1; k < p_max; ++k, tail = t)
    {
        for (arma::uword i = 0; i < dx.n_cols; ++i)
        {
            arma::uword head = heads(i);
            heads(i) = t;
            for (arma::uword j = head; j < tail; ++j, ++t)
            {
                monomials(t) = dx(i) * monomials(j);
            }
        }
    }

    return monomials;
}


}
