#pragma once

#include "nchoosek.hpp"

#include <armadillo>


namespace ifgt {


inline arma::uword get_p_max_total(arma::uword d, arma::uword p_max) {
    return nchoosek(p_max - 1 + d, d);
}
}
