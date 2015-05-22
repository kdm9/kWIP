# Libary
ADD_LIBRARY(libkmerclust
            kmerclust.cc
            countmin.cc
            kernel.cc
            population.cc
            kernels/d2.cc
            kernels/d2thresh.cc
            kernels/d2pop.cc
            kernels/d2ent.cc
            kernels/js.cc
            )

SET_TARGET_PROPERTIES(libkmerclust PROPERTIES OUTPUT_NAME kmerclust)
TARGET_LINK_LIBRARIES(libkmerclust ${KMERCLUST_DEPENDS_LIBS})

# Executable
ADD_EXECUTABLE(kmerclust kmerclust-cli.cc)
TARGET_LINK_LIBRARIES(kmerclust ${KMERCLUST_DEPENDS_LIBS} libkmerclust)
INSTALL(TARGETS kmerclust DESTINATION "bin")
