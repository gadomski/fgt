cmake_minimum_required(VERSION 2.8.12)
set(PROJECT_VERSION 0.4.5)
set(PROJECT_SOVERSION 0)

if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW) # Project version
    project(fgt LANGUAGES CXX C VERSION ${PROJECT_VERSION})
else()
    project(fgt CXX C)
endif()

message(STATUS "[fgt] Compiler: ${CMAKE_CXX_COMPILER_ID}")

configure_file(src/version.cpp.in
    ${PROJECT_BINARY_DIR}/version.cpp
    )
message(STATUS "[fgt] Version: ${PROJECT_VERSION}")

option(WITH_BENCH "Build benchmark executable" OFF)
option(WITH_OPENMP "Use OpenMP parallelization" OFF)
option(WITH_TESTS "Build test suite" ON)
message(STATUS "[fgt] With OpenMP: ${WITH_OPENMP}")

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW) # MACOSX_RPATH
endif()
if(POLICY CMP0054)
    cmake_policy(SET CMP0054 NEW) # Quotes variables in IF statments
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
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_SOVERSION}
    )
if(WITH_OPENMP)
    target_compile_definitions(Library-C++
        PUBLIC
        FGT_WITH_OPENMP
        )
endif()
if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    target_compile_options(Library-C++
        PUBLIC
        -std=c++11
        PRIVATE
        -Wall
        -pedantic
        # The no-nested-anon-types is for nanoflann
        -Wno-nested-anon-types
        )
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    target_compile_options(Library-C++
        PUBLIC
        -std=c++11
        # For some deprecations inside Eigen
        -Wno-deprecated-declarations
        PRIVATE
        -Wall
        -pedantic
        # For pragma omp's
        -Wno-unknown-pragmas
        )
elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "MSVC")
    target_compile_options(Library-C++
        PUBLIC
        # Eigen does relative includes
        /wd4464
        # This is the famous "unreferenced inline function has been removed"
        /wd4514
        # Yo, I don't care if that function wasn't inlined
        /wd4710
        # Padding
        /wd4820
        # Automatic inline expansion
        /wd4711
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

include(GenerateExportHeader)
generate_export_header(Library-C++)

include(CMakePackageConfigHelpers)
configure_file(cmake/fgt-config.cmake
    ${PROJECT_BINARY_DIR}/fgt-config.cmake
    @ONLY
    )
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/fgt-config-version.cmake
    VERSION ${PROJECT_VERSION}
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
