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

set(deps
    armadillo
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${PROJECT_SOURCE_DIR}/vendor/nanoflann-v1.1.8-25-g6c91a0e/include
    ${ARMADILLO_INCLUDE_DIRS}
    )

add_library(fgt-lib ${src})
target_link_libraries(fgt-lib ${deps})
fgt_configure_target(fgt-lib NAME fgt)
install(
    TARGETS fgt-lib
    EXPORT ${FGT_EXPORT}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
