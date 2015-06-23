# Libary
ADD_LIBRARY(libkwip
            kwip.cc
            countmin.cc
            kernel.cc
            population.cc
            kernels/ip.cc
            kernels/wip.cc
            )

SET_TARGET_PROPERTIES(libkwip PROPERTIES OUTPUT_NAME kwip)
TARGET_LINK_LIBRARIES(libkwip ${KMERCLUST_DEPENDS_LIBS})

# Executable
ADD_EXECUTABLE(kmerclust
               kmerclust-cli.cc
               kernels/js.cc
               kernels/d2freq.cc
               kernels/d2thresh.cc
               kernels/d2pop.cc
               )
TARGET_LINK_LIBRARIES(kmerclust ${KMERCLUST_DEPENDS_LIBS} libkwip)
INSTALL(TARGETS kmerclust DESTINATION "bin")

# Executable
ADD_EXECUTABLE(kwip kwip-cli.cc)
TARGET_LINK_LIBRARIES(kwip ${KMERCLUST_DEPENDS_LIBS} libkwip)
INSTALL(TARGETS kwip DESTINATION "bin")
