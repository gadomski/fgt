project(fgt CXX C)
cmake_minimum_required(VERSION 2.8.12)

set(FGT_VERSION 0.4.0)
set(FGT_SOVERSION 0)
configure_file(src/version.cpp.in
    ${PROJECT_BINARY_DIR}/version.cpp
    )

option(WITH_BENCH "Build benchmark executable" OFF)
option(WITH_OPENMP "Use OpenMP parallelization" OFF)
option(WITH_TESTS "Build test suite" ON)
message(STATUS "[fgt] With OpenMP: ${WITH_OPENMP}")

set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW) # MACOSX_RPATH
endif()

include(GetGitRevisionDescription)
git_describe(fgt_git_version --abbrev=0)
message(STATUS "[fgt] Version: ${FGT_VERSION}")
if(fgt_git_version)
    string(SUBSTRING ${fgt_git_version} 1 -1 fgt_git_version)
    if(NOT FGT_VERSION VERSION_EQUAL fgt_git_version)
        message(FATAL_ERROR "[fgt] Version (${FGT_VERSION}) does not match latest git tag (${fgt_git_version})")
    else()
        git_describe(fgt_git_describe)
        message(STATUS "[fgt] Git description: ${fgt_git_describe}")
    endif()
endif()

find_package(Eigen3 3.1 REQUIRED)
message(STATUS "[fgt] Eigen3 version: ${EIGEN3_VERSION}")

set(GTEST_SOURCE_DIR ${PROJECT_SOURCE_DIR}/vendor/gtest-1.7.0)
set(NANOFLANN_SOURCE_DIR ${PROJECT_SOURCE_DIR}/vendor/nanoflann-1.1.9)

set(cluster_src src/cluster-sequential.cpp)
if(WITH_OPENMP)
    find_package(OpenMP REQUIRED)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(cluster_src src/cluster-openmp.cpp)
endif()

add_library(Library-C++
    src/cluster.cpp
    src/direct.cpp
    src/direct_tree.cpp
    src/ifgt.cpp
    src/openmp.cpp
    src/transform.cpp
    ${cluster_src}
    ${PROJECT_BINARY_DIR}/version.cpp
    )
target_include_directories(Library-C++
    PUBLIC
    ${EIGEN3_INCLUDE_DIR}
    INTERFACE
    ${CMAKE_INSTALL_PREFIX}/include
    PRIVATE
    include
    ${NANOFLANN_SOURCE_DIR}/include
    )
set_target_properties(Library-C++ PROPERTIES
    OUTPUT_NAME fgt
    VERSION ${FGT_VERSION}
    SOVERSION ${FGT_SOVERSION}
    )
if(WITH_OPENMP)
    target_compile_definitions(Library-C++
        PUBLIC
        FGT_WITH_OPENMP
        )
endif()
target_compile_options(Library-C++
    PRIVATE
    # The no-nested-anon-types is for nanoflann
    -pedantic -Wall -Wno-nested-anon-types
    )
if(${CMAKE_COMPILER_IS_GNUCXX})
    target_compile_options(Library-C++
        PUBLIC
        # For some deprecations inside Eigen
        -Wno-deprecated-declarations
        PRIVATE
        # For pragma omp's
        -Wno-unknown-pragmas
        )
endif()
install(FILES
    include/fgt.hpp
    DESTINATION include
    )
install(TARGETS Library-C++
    DESTINATION lib
    EXPORT fgt-targets
    )

configure_file(Doxyfile.in ${PROJECT_BINARY_DIR}/Doxyfile)
add_custom_target(docs
    COMMAND doxygen ${PROJECT_BINARY_DIR}/Doxyfile
    )

if(WITH_BENCH)
    add_subdirectory(bench)
endif()

if(WITH_TESTS)
    enable_testing()
    add_subdirectory(${GTEST_SOURCE_DIR})
    set_target_properties(gtest PROPERTIES MACOSX_RPATH ON)
    set_target_properties(gtest_main PROPERTIES MACOSX_RPATH ON)

    add_subdirectory(test)
endif()

include(CMakePackageConfigHelpers)
configure_file(cmake/fgt-config.cmake.in
    ${PROJECT_BINARY_DIR}/fgt-config.cmake
    )
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/fgt-config-version.cmake
    VERSION ${FGT_VERSION}
    COMPATIBILITY AnyNewerVersion
    )
install(FILES
    ${PROJECT_BINARY_DIR}/fgt-config.cmake
    ${PROJECT_BINARY_DIR}/fgt-config-version.cmake
    DESTINATION lib/cmake/fgt
    )
install(EXPORT fgt-targets
    DESTINATION lib/cmake/fgt
    NAMESPACE Fgt::
    )
