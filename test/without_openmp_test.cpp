#include "gtest/gtest.h"

#include "fgt.hpp"

namespace fgt {

TEST(OpenMP, ApiCall) { EXPECT_FALSE(fgt::with_openmp()); }
}
