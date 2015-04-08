# Top-level CMakeLists.txt for Fgt

cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(fgt CXX C)
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/lib)

set(ARMA_MINIMUM_VERSION 4.600.1)


# Version
include(GetGitRevisionDescription)
git_describe(FGT_VERSION_DESCRIBE)
string(SUBSTRING ${FGT_VERSION_DESCRIBE} 1 -1 FGT_VERSION)
set(FGT_SOVERSION 0)

message(STATUS "fgt version ${FGT_VERSION}")


# Options
option(ARMA_64BIT_WORD "Build 64 bit Armadillo words" ON)
option(ARMA_NO_DEBUG "Build without Armadillo bounds checks" OFF)
option(BUILD_SHARED_LIBS "Create shared libraries" ON)
option(BUILD_TESTS "Build test suite" ON)
option(WITH_OPENMP "Build with OpenMP paralellization support" OFF)

message(STATUS "fgt options:
        ARMA_64BIT_WORD: ${ARMA_64BIT_WORD}
        ARMA_NO_DEBUG: ${ARMA_NO_DEBUG}
        WITH_OPENMP: ${WITH_OPENMP}")


# RPath
set(CMAKE_MACOSX_RPATH ON)
set(CMAKE_SKIP_BUILD_RPATH  OFF)
set(CMAKE_BUILD_WITH_INSTALL_RPATH OFF)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH ON)
list(FIND CMAKE_PLATFORM_IMPLICIT_LINK_DIRECTORIES "${CMAKE_INSTALL_PREFIX}/lib" isSystemDir)
if("${isSystemDir}" STREQUAL "-1")
   set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")
endif()


# Dependencies
# NO_CMAKE_SYSTEM_PATH because cmake's built-in arma find script needs to
# be superceded by the installed version
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
    endif()
    message(STATUS "OpenMP library: ${OpenMP_LIBRARY}")
endif()


# Subdirectories
include(configure_target)
add_subdirectory(include)
add_subdirectory(src)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(vendor/gtest-1.7.0)
    add_subdirectory(test)
endif()


# FgtConfig.cmake
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
