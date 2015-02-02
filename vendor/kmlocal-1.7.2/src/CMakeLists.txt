set(src
    KM_ANN.cpp
    KMeans.cpp
    KMterm.cpp
    KMrand.cpp
    KCutil.cpp
    KCtree.cpp
    KMdata.cpp
    KMcenters.cpp
    KMfilterCenters.cpp
    KMlocal.cpp
    )

include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-register")

add_library(${KMLOCAL_LIB_TARGET} STATIC ${src})
set_property(TARGET ${KMLOCAL_LIB_TARGET}
    PROPERTY OUTPUT_NAME ${KMLOCAL_LIB_NAME})
