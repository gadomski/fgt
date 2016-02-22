#include "gtest/gtest.h"

#include "fgt.hpp"

namespace fgt {

TEST(Version, Compiles) { std::string v = version(); }

TEST(GitDescribe, Compiles) { std::string d = git_describe(); }
}
