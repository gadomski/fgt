cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(fgt CXX C)

cmake_policy(SET CMP0042 NEW) # MACOSX_RPATH

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

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
    # We manually set the flags for clang-omp
    # CMake comes with a FindOpenMP script, but it doesn't work out of the box
    # on OSX with a custom-built clang-omp. Since that's what I have, I get to
    # set it the way that I want. In the future, when we use this on other OSes,
    # this will have to be smartened up.
    #
    # Note though that this works for my Travis as well.
    set(OpenMP_C_FLAGS -fopenmp CACHE STRING "OpenMP C flag")
    set(OpenMP_CXX_FLAGS -fopenmp CACHE STRING "OpenMP CXX flag")

    # This os-specific stuff is hackey, but it works for now. At the very least
    # we should use find_library for GNU as well, but building up a Travis simulation
    # box with openmp is beyond my scope of work ATM :-).
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(OpenMP_LIBRARY gomp CACHE STRING "OpenMP library name (usually gomp for gcc)")
    else()
        find_library(OpenMP_LIBRARY iomp5)
        find_path(OpenMP_INCLUDE_DIR omp.h
            PATH_SUFFIXES libiomp
            )
        message(STATUS "OpenMP include directory: ${OpenMP_INCLUDE_DIR}")
    endif()
    message(STATUS "OpenMP library: ${OpenMP_LIBRARY}")

endif()

include(configure_target)
add_subdirectory(include)
add_subdirectory(src)

if(WITH_TESTS)
    enable_testing()
    add_subdirectory(test)

    add_subdirectory(vendor/gtest-1.7.0)
    set_target_properties(gtest PROPERTIES
        MACOSX_RPATH ON
        )
    set_target_properties(gtest_main PROPERTIES
        MACOSX_RPATH ON
        )
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
