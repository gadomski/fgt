add_executable(bench bench.cpp)
target_include_directories(bench
    PRIVATE
        ${PROJECT_SOURCE_DIR}/include
    )
target_link_libraries(bench
    PRIVATE
        Library-C++
    )
