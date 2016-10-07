INCLUDE_DIRECTORIES(ext/greatest/
                    ${CMAKE_CURRENT_SOURCE_DIR})

ADD_EXECUTABLE(kwip_tests ${CMAKE_CURRENT_SOURCE_DIR}/main.c)
TARGET_LINK_LIBRARIES(kwip_tests ${KWIP_DEP_LIBS} kwiplib blosc_static)
