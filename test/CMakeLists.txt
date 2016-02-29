function(fgt_test target)
    add_executable(${target} ${ARGN})
    target_include_directories(${target}
        PRIVATE
        ${GTEST_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
        ${PROJECT_BINARY_DIR}
        )
    target_link_libraries(${target}
        PRIVATE
        gtest_main
        Library-C++
        )
    if(${CMAKE_COMPILER_IS_GNUCXX})
        target_compile_options(${target}
            PRIVATE
            # For pragma omp's
            -Wno-unknown-pragmas
            )
    endif()
    add_test(NAME ${target} COMMAND ${target})
endfunction()

configure_file(
    support.hpp.in
    ${CMAKE_CURRENT_BINARY_DIR}/support.hpp
    )

fgt_test(cluster cluster_test.cpp)
fgt_test(direct direct_test.cpp)
fgt_test(direct-tree direct_tree_test.cpp)
fgt_test(ifgt ifgt_test.cpp)
fgt_test(version version_test.cpp)

if(WITH_OPENMP)
    fgt_test(with_openmp with_openmp_test.cpp)
else()
    fgt_test(without_openmp without_openmp_test.cpp)
endif()
