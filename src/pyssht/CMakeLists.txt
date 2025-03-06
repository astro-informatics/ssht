add_custom_command(
    OUTPUT cpyssht.c
    MAIN_DEPENDENCY cpyssht.pyx
    VERBATIM
    COMMAND Python::Interpreter -m cython "${CMAKE_CURRENT_SOURCE_DIR}/cpyssht.pyx"
            --output-file "${CMAKE_CURRENT_BINARY_DIR}/cpyssht.c")
add_custom_command(
    OUTPUT ducc_interface.c
    MAIN_DEPENDENCY ducc_interface.pyx
    VERBATIM
    COMMAND
        Python::Interpreter -m cython "${CMAKE_CURRENT_SOURCE_DIR}/ducc_interface.pyx"
        --output-file "${CMAKE_CURRENT_BINARY_DIR}/ducc_interface.c")

python_add_library(cpyssht MODULE "${CMAKE_CURRENT_BINARY_DIR}/cpyssht.c" WITH_SOABI)
python_add_library(ducc_interface MODULE "${CMAKE_CURRENT_BINARY_DIR}/ducc_interface.c"
                   WITH_SOABI)
target_link_libraries(cpyssht PUBLIC ssht)
target_include_directories(cpyssht PUBLIC ${Python3_NumPy_INCLUDE_DIR})
target_include_directories(ducc_interface PUBLIC ${Python3_NumPy_INCLUDE_DIR})
if(WIN32)
    target_compile_definitions(cpyssht PUBLIC MS_WIN64)
    target_compile_definitions(ducc_interface PUBLIC MS_WIN64)
endif()
install(TARGETS cpyssht LIBRARY DESTINATION ${SKBUILD_PROJECT_NAME})
