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


#ADD_EXECUTABLE(kmerclust
#               kmerclust/kmerclust-cli.cc
#               kmerclust/js.cc
#               kmerclust/d2freq.cc
#               kmerclust/d2thresh.cc
#               kmerclust/d2pop.cc
#               )
#TARGET_LINK_LIBRARIES(kmerclust ${KMERCLUST_DEPENDS_LIBS} libkwip)


ADD_EXECUTABLE(kwip kwip-cli.cc)
TARGET_LINK_LIBRARIES(kwip ${KMERCLUST_DEPENDS_LIBS} libkwip)
INSTALL(TARGETS kwip DESTINATION "bin")


ADD_EXECUTABLE(kwip-entvec utils/entvec.cc)
TARGET_LINK_LIBRARIES(kwip-entvec ${KMERCLUST_DEPENDS_LIBS} libkwip)
INSTALL(TARGETS kwip-entvec DESTINATION "bin")
