set(src
    clustering.cpp
    constant_series.cpp
    direct.cpp
    direct_tree.cpp
    fgt.cpp
    gauss_transform.cpp
    ifgt.cpp
    monomials.cpp
    truncation_number.cpp
    )

add_library(Library ${src})
configure_target(Library)
set_target_properties(Library PROPERTIES
    VERSION ${FGT_VERSION}
    SOVERSION ${FGT_SOVERSION}
    OUTPUT_NAME fgt
    )

target_include_directories(Library
    PUBLIC
    ${ARMADILLO_INCLUDE_DIRS}
    PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/vendor/nanoflann-v1.1.8-25-g6c91a0e/include
    )

target_link_libraries(Library PUBLIC armadillo)
if (WITH_OPENMP)
    target_compile_definitions(Library PRIVATE FGT_WITH_OPENMP)
endif()

install(
    TARGETS Library
    EXPORT FgtTargets
    INCLUDES DESTINATION include
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
