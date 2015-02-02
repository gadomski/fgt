set(src
    choose_parameters.cpp
    choose_truncation_number.cpp
    ifgt.cpp
    )

set(deps
    ${ARMADILLO_LIBRARIES}
    )

include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${ARMADILLO_INCLUDE_DIRS}
    )

add_library(${IFGT_LIB_TARGET} ${src})
target_link_libraries(${IFGT_LIB_TARGET} ${deps})
set_property(TARGET ${IFGT_LIB_TARGET}
    PROPERTY OUTPUT_NAME ${IFGT_LIB_NAME})
