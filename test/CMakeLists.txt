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

add_executable(Test ${src})
configure_target(Test)
set_target_properties(Test PROPERTIES OUTPUT_NAME fgt-test)
target_include_directories(Test
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_BINARY_DIR}
        ${gtest_SOURCE_DIR}/include
    )
target_link_libraries(Test PRIVATE Library-C++ gtest_main)
add_test(NAME unit-tests COMMAND Test)

add_executable(Trial trial.cpp)
configure_target(Trial)
set_target_properties(Trial PROPERTIES OUTPUT_NAME fgt-trial)
target_include_directories(Trial
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
        ${gtest_SOURCE_DIR}/include
    )
target_link_libraries(Trial PRIVATE Library-C++ gtest_main)
