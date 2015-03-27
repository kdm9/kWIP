# Libary
ADD_LIBRARY(libkmerclust distance.cc
                         d2.cc
                         d2pop.cc
                         population.cc
                         countmin.cc)
SET_TARGET_PROPERTIES(libkmerclust PROPERTIES OUTPUT_NAME kmerclust)

# Executable
ADD_EXECUTABLE(klust main.cc)
TARGET_LINK_LIBRARIES(klust ${DEPENDS_LIBS} libkmerclust)
INSTALL(TARGETS klust DESTINATION "bin")
