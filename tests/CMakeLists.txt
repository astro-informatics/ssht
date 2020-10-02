add_executable(test_ssht test_ssht.c utilities.c interface.c)
target_link_libraries(test_ssht ssht cmocka)
set_target_properties(test_ssht PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                           ${PROJECT_BINARY_DIR} C_STANDARD 11)
add_test(NAME test_ssht COMMAND test_ssht)
