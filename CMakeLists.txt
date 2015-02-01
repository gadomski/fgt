# Basic setup
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(ifgt CXX C)


# Options
option(BUILD_TESTS "Build test suite" ON)


# Targets
set(IFGT_TEST_TARGET ifgt-test)


# Subdirectories
if(BUILD_TESTS)
    add_subdirectory(test)
    add_subdirectory(vendor/gtest-1.7.0)
endif()
