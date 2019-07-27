SET(DEACTIVATE_LZ4 ON CACHE BOOL "")
SET(DEACTIVATE_SNAPPY ON CACHE BOOL "")
SET(DEACTIVATE_ZSTD OFF CACHE BOOL "")
SET(DEACTIVATE_ZLIB ON CACHE BOOL "")
SET(PREFER_EXTERNAL_ZSTD OFF CACHE BOOL "")
SET(BUILD_SHARED OFF CACHE BOOL "")
SET(BUILD_TESTS OFF CACHE BOOL "")
SET(BUILD_BENCHMARKS OFF CACHE BOOL "")
ADD_SUBDIRECTORY(ext/blosc)
set_target_properties(blosc_static PROPERTIES COMPILE_FLAGS -w)

FILE(GLOB DEP_SRCS
    ext/blosc-h5/src/blosc_*.c
    ext/xxhash.c
    ext/clogged/clogged.c
)
set_source_files_properties(${DEP_SRCS} PROPERTIES COMPILE_FLAGS -w)

INCLUDE_DIRECTORIES(SYSTEM
    ext/blosc/blosc
    ext/blosc-h5/src
    ext/
)


SET(KWIP_SRCS
    kwip_array.c
    kwip_distcalc.c
    kwip_kmercount.c
    kwip_kmerhash.c
    kwip_metrics.c
    kwip_omp_dist.c
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
