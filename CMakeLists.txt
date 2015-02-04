# Basic setup
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(ifgt CXX C)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(CMAKE_MACOSX_RPATH ON)

set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/lib)


# Options
option(BUILD_TESTS "Build test suite" ON)


# Dependencies
find_package(Armadillo 4.500.1 REQUIRED)


# Targets
set(IFGT_LIB_TARGET ifgt-lib)
set(IFGT_LIB_NAME ifgt)
set(IFGT_TEST_TARGET ifgt-test)


# Subdirectories
add_subdirectory(include)
add_subdirectory(src)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(vendor/gtest-1.7.0)
    add_subdirectory(test)
endif()
