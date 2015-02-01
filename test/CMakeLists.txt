set(src
    choose_parameters_test.cpp
    )

set(deps
    gtest_main)


add_executable(${IFGT_TEST_TARGET} ${src})
target_link_libraries(${IFGT_TEST_TARGET} ${deps})
