set(src
    choose_parameters_test.cpp
    choose_truncation_number_test.cpp
    clusters_test.cpp
    )

set(deps
    ${IFGT_LIB_TARGET}
    gtest_main
    ${KMLOCAL_LIB_TARGET}
    )

include_directories(
    ${gtest_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${kmlocal_SOURCE_DIR}/src
    )

add_executable(${IFGT_TEST_TARGET} ${src})
target_link_libraries(${IFGT_TEST_TARGET} ${deps})
add_test(unit-tests ${PROJECT_BINARY_DIR}/bin/${IFGT_TEST_TARGET})
