# Basic setup
cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(ifgt CXX C)
set(IFGT_VERSION 0.1.0)

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
set(IFGT_EXPORT IfgtTargets)


# Subdirectories
add_subdirectory(include)
add_subdirectory(src)

if(BUILD_TESTS)
    enable_testing()
    add_subdirectory(vendor/gtest-1.7.0)
    add_subdirectory(test)
endif()


# IfgtConfig.cmake
include(CMakePackageConfigHelpers)

set(INCLUDE_INSTALL_DIR include)

configure_package_config_file(
    cmake/IfgtConfig.cmake.in
    ${PROJECT_BINARY_DIR}/IfgtConfig.cmake
    INSTALL_DESTINATION lib/ifgt/cmake
    PATH_VARS INCLUDE_INSTALL_DIR
    )
write_basic_package_version_file(
    ${PROJECT_BINARY_DIR}/IfgtConfigVersion.cmake
    VERSION ${IFGT_VERSION}
    COMPATIBILITY SameMajorVersion
    )
export(TARGETS ${IFGT_LIB_TARGET} FILE ${PROJECT_BINARY_DIR}/IfgtTargets.cmake)

install(
    FILES
        ${PROJECT_BINARY_DIR}/IfgtConfig.cmake
        ${PROJECT_BINARY_DIR}/IfgtConfigVersion.cmake
    DESTINATION
        lib/ifgt/cmake
        )
install(
    EXPORT ${IFGT_EXPORT}
    DESTINATION lib/ifgt/cmake
    )
