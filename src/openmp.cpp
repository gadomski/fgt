#include "fgt.hpp"

namespace fgt {

bool with_openmp() {
#ifdef FGT_WITH_OPENMP
    return true;
#else
    return false;
#endif
}
}
