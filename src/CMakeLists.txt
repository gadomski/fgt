set(src
    choose_parameters.cpp
    choose_truncation_number.cpp
    )


add_library(${IFGT_LIB_TARGET} ${src})
set_property(TARGET ${IFGT_LIB_TARGET}
    PROPERTY OUTPUT_NAME ${IFGT_LIB_NAME})
