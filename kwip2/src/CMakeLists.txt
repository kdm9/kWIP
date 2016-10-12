SET(DEACTIVATE_SNAPPY ON CACHE BOOL "")
SET(DEACTIVATE_ZSTD ON CACHE BOOL "")
SET(DEACTIVATE_ZLIB ON CACHE BOOL "")
SET(BUILD_SHARED OFF CACHE BOOL "")
SET(BUILD_TESTS OFF CACHE BOOL "")
SET(BUILD_BENCHMARKS OFF CACHE BOOL "")
ADD_SUBDIRECTORY(ext/blosc)
INCLUDE_DIRECTORIES(ext/blosc/blosc)

FILE(GLOB DEP_SRCS
     ext/blosc-h5/src/blosc_*.c
     ext/xxhash.c
     ext/clogged/clogged.c
     ext/commander/commander.c
     )

INCLUDE_DIRECTORIES(ext/blosc-h5/src/
                    ext/commander/
                    ext/clogged/
                    ext/
                    )


SET(KWIP_SRCS
    kwip_array.c
    kwip_kernelcalc.c
    kwip_kmercount.c
    kwip_kmerhash.c
    kwip_utils.c
)

ADD_LIBRARY(kwiplib STATIC ${KWIP_SRCS} ${DEP_SRCS})
TARGET_LINK_LIBRARIES(kwiplib blosc_static ${KWIP_DEP_LIBS})
SET_TARGET_PROPERTIES(kwiplib PROPERTIES OUTPUT_NAME kwip)

ADD_EXECUTABLE(kwip omp_main.c)
TARGET_LINK_LIBRARIES(kwip kwiplib ${KWIP_DEP_LIBS} m)

add_subdirectory(tests)
