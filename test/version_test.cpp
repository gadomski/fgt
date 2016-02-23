#include "gtest/gtest.h"

#include "test/support.hpp"

namespace fgt {

TEST(Version, Compiles) { std::string v = version(); }
}
