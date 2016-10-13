
ADD_EXECUTABLE(kwip_tests ${CMAKE_CURRENT_SOURCE_DIR}/main.c)
TARGET_LINK_LIBRARIES(kwip_tests kwiplib ${CMOCKA_LIBRARIES})
TARGET_INCLUDE_DIRECTORIES(kwip_tests SYSTEM PRIVATE ${CMOCKA_INCLUDE_DIRS})

ADD_TEST(NAME run_tests
         COMMAND ${CMAKE_BINARY_DIR}/bin/kwip_tests
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

# Copy test files over to bin dir
ADD_CUSTOM_COMMAND(TARGET kwip_tests
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
		                ${CMAKE_CURRENT_SOURCE_DIR}/data
                        ${CMAKE_BINARY_DIR}/data)
