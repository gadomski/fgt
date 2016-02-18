cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(fgt CXX C)

if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW) # MACOSX_RPATH
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(ConfigureTarget)
include(GetGitRevisionDescription)

git_describe(FGT_VERSION_DESCRIBE)
if(NOT FGT_VERSION_DESCRIBE STREQUAL "HEAD-HASH-NOTFOUND")
    string(SUBSTRING ${FGT_VERSION_DESCRIBE} 1 -1 FGT_VERSION)
else()
    file(STRINGS VERSION.txt FGT_VERSION)
    string(SUBSTRING ${FGT_VERSION} 1 -1 FGT_VERSION)
endif()
set(FGT_SOVERSION 0)
message(STATUS "fgt version ${FGT_VERSION}")

option(ARMA_64BIT_WORD "Build 64 bit Armadillo words" ON)
option(ARMA_NO_DEBUG "Build without Armadillo bounds checks" OFF)
option(WITH_TESTS "Build test suite" ON)
option(WITH_OPENMP "Build with OpenMP paralellization support" OFF)

message(STATUS "fgt options:
    ARMA_64BIT_WORD: ${ARMA_64BIT_WORD}
      ARMA_NO_DEBUG: ${ARMA_NO_DEBUG}
         WITH_TESTS: ${WITH_TESTS}
        WITH_OPENMP: ${WITH_OPENMP}")

set(ARMA_MINIMUM_VERSION 4.600.1)
# We set NO_CMAKE_SYSTEM_PATH because newer Armadillos come with a configure
# script that conflicts with cmake's pre-packaged version
find_package(Armadillo ${ARMA_MINIMUM_VERSION} REQUIRED NO_CMAKE_SYSTEM_PATH)

if (WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_library(Library-C++
    src/clustering.cpp
    src/constant_series.cpp
    src/direct.cpp
    src/direct_tree.cpp
    src/fgt.cpp
    src/gauss_transform.cpp
    src/ifgt.cpp
    src/monomials.cpp
    src/truncation_number.cpp
    )
configure_target(Library-C++)
target_include_directories(Library-C++
    PUBLIC
    ${ARMADILLO_INCLUDE_DIRS}
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/vendor/nanoflann-v1.1.8-25-g6c91a0e/include
    )
target_link_libraries(Library-C++ PUBLIC armadillo)
set_target_properties(Library-C++ PROPERTIES
    VERSION ${FGT_VERSION}
    SOVERSION ${FGT_SOVERSION}
    OUTPUT_NAME fgt
    )
if (WITH_OPENMP)
    target_compile_definitions(Library-C++ PRIVATE FGT_WITH_OPENMP)
endif()
install(
    TARGETS Library-C++
    EXPORT FgtTargets
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )

install(
    DIRECTORY include/fgt
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
    )

if(WITH_TESTS)
    enable_testing()
    add_subdirectory(vendor/gtest-1.7.0)
    set_target_properties(gtest PROPERTIES
        MACOSX_RPATH ON
        )
    set_target_properties(gtest_main PROPERTIES
        MACOSX_RPATH ON
        )

    add_subdirectory(test)
endif()

include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/FgtConfig.cmake.in
    ${PROJECT_BINARY_DIR}/FgtConfig.cmake
    INSTALL_DESTINATION lib/fgt/cmake
    NO_SET_AND_CHECK_MACRO
    NO_CHECK_REQUIRED_COMPONENTS_MACRO
    )
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/FgtConfigVersion.cmake
    VERSION ${FGT_VERSION}
    COMPATIBILITY SameMajorVersion
    )
install(
    FILES
        ${PROJECT_BINARY_DIR}/FgtConfig.cmake
        ${PROJECT_BINARY_DIR}/FgtConfigVersion.cmake
    DESTINATION lib/fgt/cmake
    )
install(
    EXPORT FgtTargets
    NAMESPACE Fgt::
    DESTINATION lib/fgt/cmake
    )
