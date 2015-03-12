# Basic setup
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(fgt CXX C)

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/lib)
set(FGT_EXPORT FgtTargets)

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
find_package(Armadillo 4.650.1 REQUIRED NO_CMAKE_SYSTEM_PATH)

if (WITH_OPENMP)
    find_package(OpenMP REQUIRED QUIET)
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

set(INCLUDE_INSTALL_DIR include)

configure_package_config_file(
    cmake/FgtConfig.cmake.in
    ${PROJECT_BINARY_DIR}/FgtConfig.cmake
    INSTALL_DESTINATION lib/fgt/cmake
    PATH_VARS INCLUDE_INSTALL_DIR
    )
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/FgtConfigVersion.cmake
    VERSION ${FGT_VERSION}
    COMPATIBILITY SameMajorVersion
    )
export(TARGETS fgt-lib FILE ${PROJECT_BINARY_DIR}/FgtTargets.cmake)

install(
    FILES
        ${PROJECT_BINARY_DIR}/FgtConfig.cmake
        ${PROJECT_BINARY_DIR}/FgtConfigVersion.cmake
    DESTINATION
        lib/fgt/cmake
        )
install(
    EXPORT ${FGT_EXPORT}
    DESTINATION lib/fgt/cmake
    )
