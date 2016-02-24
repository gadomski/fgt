#include "gtest/gtest.h"

#include "fgt.hpp"

namespace fgt {

TEST(OpenMP, ApiCall) { EXPECT_TRUE(fgt::with_openmp()); }
}
