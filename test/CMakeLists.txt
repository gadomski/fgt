configure_file(
    config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.hpp
    )

set(src
    clustering_test.cpp
    constant_series_test.cpp
    direct_test.cpp
    direct_tree_test.cpp
    fgt_test.cpp
    gauss_transform_test.cpp
    ifgt_test.cpp
    monomials_test.cpp
    truncation_number_test.cpp
    )

set(deps
    fgt-lib
    gtest_main
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}
    ${gtest_SOURCE_DIR}/include
    ${ARMADILLO_INCLUDE_DIRS}
    )

add_executable(fgt-test ${src})
target_link_libraries(fgt-test ${deps})
fgt_configure_target(fgt-test)
add_test(unit-tests ${PROJECT_BINARY_DIR}/bin/fgt-test)

add_executable(fgt-trial trial.cpp)
target_link_libraries(fgt-trial ${deps})
fgt_configure_target(fgt-trial)
