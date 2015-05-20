CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
PROJECT(kmerclust)

###############
# Cmake Blerg #
###############

# Find the modules
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# Set sane output paths
SET(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/lib)


ENABLE_TESTING()

# Default to Release builds
IF (NOT CMAKE_BUILD_TYPE)
	SET(CMAKE_BUILD_TYPE Release)
ENDIF()

# Version jiggery-pokery
IF (EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/version")
	# Use version string from ./version
        FILE(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/version" VERSION)
ELSE()
        # git describe as versioning
        EXECUTE_PROCESS(COMMAND git describe --always
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        OUTPUT_VARIABLE VERSION
                        OUTPUT_STRIP_TRAILING_WHITESPACE)
ENDIF()

MESSAGE(STATUS "${CMAKE_BUILD_TYPE} build of ${PROJECT_NAME} version: ${VERSION}")

############################
## User-definable options ##
############################

OPTION(ENABLE_ASAN "Build with AddressSanitizer" OFF)
OPTION(UNSAFE_MATH "Use -ffast-math in Release build? (costs a little accuracy)" ON)
OPTION(ENABLE_TSAN "Build with ThreadSanitizer" OFF)
OPTION(USE_OPENMP "Use OpenMP for parallelism" ON)

###############################
## Find Packages and Headers ##
###############################

FIND_PACKAGE(khmer 1.3 REQUIRED)
SET(KMERCLUST_DEPENDS_LIBS ${KMERCLUST_DEPENDS_LIBS} ${KHMER_LIBRARIES})
INCLUDE_DIRECTORIES(${KHMER_INCLUDE_DIRS})

IF(NOT KHMER_FOUND)
	MESSAGE(FATAL_ERROR "Khmer's C++ libary must be installed. "
	                    "See github.com/kdmurray91/kmerclust for instructions")
ENDIF()

IF(USE_OPENMP)
  FIND_PACKAGE(OpenMP)
  IF(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  ENDIF()
ENDIF()

##########################
## Set Compiler Options ##
##########################

# Set CFLAGS
SET(WARN_FLAGS "${WARN_FLAGS} -Wall")
SET(OPT_FLAGS "${OPT_FLAGS} -march=native")

SET(REL_OPT_FLAGS "-O3")
IF(UNSAFE_MATH)
	SET(REL_OPT_FLAGS "-ffast-math -Ofast")
ENDIF()

IF(ENABLE_ASAN AND ENABLE_TSAN)
	MESSAGE(FATAL_ERROR "Cannot enable both ThreadSanitizer and AddressSanitizer")
ENDIF()

IF(ENABLE_ASAN)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")
ENDIF()

IF(ENABLE_TSAN)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=thread -fPIC -pie")
ENDIF()

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c11 ${WARN_FLAGS} ${OPT_FLAGS}")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 ${WARN_FLAGS} ${OPT_FLAGS}")
SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${REL_OPT_FLAGS}")

# Headers occur in these dirs
INCLUDE_DIRECTORIES(${CMAKE_BINARY_DIR}
                    ${CMAKE_SOURCE_DIR}/src
                    ${CMAKE_SOURCE_DIR}/tests
                    ${CMAKE_SOURCE_DIR}/tests/ext)

# Finally, configure the config header
CONFIGURE_FILE(${CMAKE_SOURCE_DIR}/src/kmerclust-config.hh.in
               ${CMAKE_BINARY_DIR}/kmerclust-config.hh)

# Libary sources
ADD_SUBDIRECTORY(src)
ADD_SUBDIRECTORY(tests)
