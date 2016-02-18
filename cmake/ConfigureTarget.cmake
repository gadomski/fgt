function(configure_target target)
    target_compile_options(${target}
        PUBLIC -std=c++11
        PRIVATE -Wall -Werror -pedantic
        )
    if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${target}
            PRIVATE -Wno-nested-anon-types
            )
    endif()
    if(ARMA_64BIT_WORD)
        target_compile_definitions(${target} PUBLIC ARMA_64BIT_WORD)
    endif()
    if(ARMA_NO_DEBUG)
        target_compile_definitions(${target} PUBLIC ARMA_NO_DEBUG)
    endif()
endfunction()
