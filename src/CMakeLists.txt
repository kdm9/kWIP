# Libary
ADD_LIBRARY(libkmerclust
            countmin.cc
            distance.cc
            d2.cc
            d2thresh.cc
            d2pop.cc
            population.cc
            js.cc
            )

SET_TARGET_PROPERTIES(libkmerclust PROPERTIES OUTPUT_NAME kmerclust)
TARGET_LINK_LIBRARIES(libkmerclust ${DEPENDS_LIBS})

# Executable
ADD_EXECUTABLE(kmerclust kmerclust-cli.cc)
TARGET_LINK_LIBRARIES(kmerclust ${DEPENDS_LIBS} libkmerclust)
INSTALL(TARGETS kmerclust DESTINATION "bin")
