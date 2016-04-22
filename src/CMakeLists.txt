IF (NOT USE_SYSTEM_KHMER)
    FILE(GLOB KHMER_SRC ext/oxli/*.cc)
ENDIF()

ADD_LIBRARY(libkwip
            kwip-utils.cc
            countmin.cc
            kernel.cc
            population.cc
            kernels/ip.cc
            kernels/wip.cc
            ${KHMER_SRC}
            )

SET_TARGET_PROPERTIES(libkwip PROPERTIES OUTPUT_NAME kwip)
TARGET_LINK_LIBRARIES(libkwip ${KMERCLUST_DEPENDS_LIBS})


ADD_EXECUTABLE(kwip kwip-cli.cc)
TARGET_LINK_LIBRARIES(kwip ${KMERCLUST_DEPENDS_LIBS} libkwip)
INSTALL(TARGETS kwip DESTINATION "bin")


ADD_EXECUTABLE(kwip-stats utils/kwip-stats.cc)
TARGET_LINK_LIBRARIES(kwip-stats ${KMERCLUST_DEPENDS_LIBS} libkwip)
INSTALL(TARGETS kwip-stats DESTINATION "bin")
