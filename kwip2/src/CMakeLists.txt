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


FILE(GLOB KWIP_SRCS ${CMAKE_CURRENT_SOURCE_DIR}/kwip_*.c)

ADD_LIBRARY(kwiplib STATIC ${KWIP_SRCS} ${DEP_SRCS})
SET_TARGET_PROPERTIES(kwiplib PROPERTIES OUTPUT_NAME kwip)

ADD_EXECUTABLE(kwip ${CMAKE_CURRENT_SOURCE_DIR}/main.c)
TARGET_LINK_LIBRARIES(kwip kwiplib blosc_static ${KWIP_DEP_LIBS})

add_subdirectory(tests)
