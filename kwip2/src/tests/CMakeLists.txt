INCLUDE_DIRECTORIES(ext/greatest/
                    ${CMAKE_CURRENT_SOURCE_DIR})

ADD_EXECUTABLE(kwip_tests ${CMAKE_CURRENT_SOURCE_DIR}/main.c)
TARGET_LINK_LIBRARIES(kwip_tests ${KWIP_DEP_LIBS} kwiplib blosc_static)

ADD_TEST(NAME run_tests
         COMMAND ${CMAKE_BINARY_DIR}/bin/kwip_tests
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Copy test files over to bin dir
ADD_CUSTOM_COMMAND(TARGET kwip_tests
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
		                ${CMAKE_CURRENT_SOURCE_DIR}/data
                        ${CMAKE_BINARY_DIR}/data)
