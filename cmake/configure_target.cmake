function(fgt_configure_target target)
    set(one_value_args NAME EXTRA_COMPILE_FLAGS)
    cmake_parse_arguments(TARGET "" "${one_value_args}" "" ${ARGN})
    set_property(TARGET ${target} APPEND_STRING PROPERTY
        COMPILE_FLAGS " -std=c++11 -Wall -pedantic -Werror -Wno-nested-anon-types ${OpenMP_CXX_FLAGS} ${TARGET_EXTRA_COMPILE_FLAGS}")

    set_target_properties(${target} PROPERTIES
        VERSION ${FGT_VERSION}
        SOVERSION ${FGT_SOVERSION}
        )

    if (TARGET_NAME)
        set_property(TARGET ${target} PROPERTY OUTPUT_NAME ${TARGET_NAME})
    endif()

    if (WITH_OPENMP)
        target_link_libraries(${target} ${OpenMP_LIBRARY})
        target_compile_definitions(${target} PRIVATE FGT_WITH_OPENMP)
    endif()
    if(ARMA_64BIT_WORD)
        target_compile_definitions(${target} PUBLIC ARMA_64BIT_WORD)
    endif()
    if(ARMA_NO_DEBUG)
        target_compile_definitions(${target} PUBLIC ARMA_NO_DEBUG)
    endif()
endfunction()
