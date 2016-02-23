#include <fstream>

#include "gtest/gtest.h"

#include "test/support.hpp"

namespace fgt {

TEST(Version, MatchesVersionFile) {
    std::string v = version();
    std::ifstream file(project_source_filename("VERSION.txt"));
    std::string version_txt;
    std::getline(file, version_txt);
    EXPECT_EQ(version_txt, v);
}
}
