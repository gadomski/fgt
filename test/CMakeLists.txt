configure_file(
    config.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/config.hpp
    )

set(src
    choose_parameters_test.cpp
    choose_truncation_number_test.cpp
    clusters_test.cpp
    )

set(deps
    ${IFGT_LIB_TARGET}
    gtest_main
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/src
    ${CMAKE_CURRENT_BINARY_DIR}
    ${gtest_SOURCE_DIR}/include
    ${kmlocal_SOURCE_DIR}/src
    )

add_executable(${IFGT_TEST_TARGET} ${src})
target_link_libraries(${IFGT_TEST_TARGET} ${deps})
add_test(unit-tests ${PROJECT_BINARY_DIR}/bin/${IFGT_TEST_TARGET})
