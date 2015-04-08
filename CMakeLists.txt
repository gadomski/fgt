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
    find_package(OpenMP REQUIRED)
    message(STATUS "OpenMP compiler detected: ${OpenMP_CXX_FLAGS}")

    # We punt on finding non-gomp OpenMP because clang is super-variable. This
    # is becuase Clang (as of this writing) doesn't support OpenMP out of the
    # box, so clang+openmp setups are probably using
    # https://clang-omp.github.io/. This will change eventually, and hopefully
    # this will all converge -- eventually.
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        set(OpenMP_LIBRARY gomp CACHE STRING "OpenMP library name (usually gomp for gcc)")
    else()
        set(OpenMP_LIBRARY OpenMP_LIBRARY-NOTFOUND CACHE STRING "OpenMP library")
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
