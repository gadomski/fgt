set(src
    clustering.cpp
    constant_series.cpp
    direct.cpp
    gaussian_transform.cpp
    ifgt.cpp
    monomials.cpp
    )

set(deps
    ${ARMADILLO_LIBRARIES}
    )

set(public_compile_defs
    ARMA_64BIT_WORD
    )

if (ARMA_NO_DEBUG)
    list(APPEND public_compile_defs ARMA_NO_DEBUG)
endif()

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${ARMADILLO_INCLUDE_DIRS}
    ${kmlocal_SOURCE_DIR}/src
    )

add_library(${FGT_LIB_TARGET} ${src})
target_link_libraries(${FGT_LIB_TARGET} ${deps})
set_property(TARGET ${FGT_LIB_TARGET}
    PROPERTY OUTPUT_NAME ${FGT_LIB_NAME})
target_compile_definitions(${FGT_LIB_TARGET} PUBLIC ${public_compile_defs})

install(
    TARGETS ${FGT_LIB_TARGET}
    EXPORT ${FGT_EXPORT}
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
    )
