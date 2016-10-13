SET(DEACTIVATE_SNAPPY ON CACHE BOOL "")
SET(DEACTIVATE_ZSTD ON CACHE BOOL "")
SET(DEACTIVATE_ZLIB ON CACHE BOOL "")
SET(BUILD_SHARED OFF CACHE BOOL "")
SET(BUILD_TESTS OFF CACHE BOOL "")
SET(BUILD_BENCHMARKS OFF CACHE BOOL "")
ADD_SUBDIRECTORY(ext/blosc)

FILE(GLOB DEP_SRCS
    ext/blosc-h5/src/blosc_*.c
    ext/xxhash.c
    ext/clogged/clogged.c
)

INCLUDE_DIRECTORIES(SYSTEM
    ext/blosc/blosc
    ext/blosc-h5/src
    ext/
)


SET(KWIP_SRCS
    kwip_array.c
    kwip_kernelcalc.c
    kwip_kmercount.c
    kwip_kmerhash.c
    kwip_metric_ip.c
    kwip_omp_kernel.c
    kwip_utils.c
)

ADD_LIBRARY(kwiplib STATIC ${KWIP_SRCS} ${DEP_SRCS})
TARGET_LINK_LIBRARIES(kwiplib blosc_static ${KWIP_DEP_LIBS})
SET_TARGET_PROPERTIES(kwiplib PROPERTIES OUTPUT_NAME kwip)

ADD_EXECUTABLE(kwip omp_main.c)
TARGET_LINK_LIBRARIES(kwip kwiplib ${KWIP_DEP_LIBS} m)

pkg_search_module(CMOCKA cmocka)

if (CMOCKA_FOUND)
    add_subdirectory(tests)
else()
    message(STATUS "CMocka not found, disabling tests")
endif()
