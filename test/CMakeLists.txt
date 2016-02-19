function(fgt_test target)
    add_executable(${target} ${ARGN})
    configure_target(${target})
    target_include_directories(${target}
        PRIVATE
            ${PROJECT_SOURCE_DIR}/include
            ${PROJECT_SOURCE_DIR}/src
            ${CMAKE_CURRENT_BINARY_DIR}
            ${gtest_SOURCE_DIR}/include
        )
    target_link_libraries(${target}
        PRIVATE
            gtest_main
            Library-C++
        )
    add_test(NAME ${target} COMMAND ${target})
endfunction()

configure_file(
    config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.hpp
    )

fgt_test(clustering clustering_test.cpp)
fgt_test(constant-series constant_series_test.cpp)
fgt_test(direct direct_test.cpp)
fgt_test(direct-tree direct_tree_test.cpp)
fgt_test(fgt fgt_test.cpp)
fgt_test(gauss-transform gauss_transform_test.cpp)
fgt_test(ifgt ifgt_test.cpp)
fgt_test(monomials monomials_test.cpp)
fgt_test(trial trial.cpp)
fgt_test(truncation-number truncation_number_test.cpp)
