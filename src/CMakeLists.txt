set(src
    choose_truncation_number.cpp
    clustering.cpp
    clustering/gonzalez.cpp
    clustering_factory.cpp
    constant_series.cpp
    ifgt.cpp
    monomials.cpp
    parameters.cpp
    )

set(deps
    ${ARMADILLO_LIBRARIES}
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${ARMADILLO_INCLUDE_DIRS}
    ${kmlocal_SOURCE_DIR}/src
    )

add_library(${IFGT_LIB_TARGET} ${src})
target_link_libraries(${IFGT_LIB_TARGET} ${deps})
set_property(TARGET ${IFGT_LIB_TARGET}
    PROPERTY OUTPUT_NAME ${IFGT_LIB_NAME})
target_compile_definitions(${IFGT_LIB_TARGET} PUBLIC ARMA_64BIT_WORD)

install(
    TARGETS ${IFGT_LIB_TARGET}
    EXPORT ${IFGT_EXPORT}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
