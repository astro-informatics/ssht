find_package(Python3 REQUIRED COMPONENTS Development Interpreter)

# look for cython
execute_process(COMMAND ${Python3_EXECUTABLE} -m cython -V
    RESULT_VARIABLE CYTHON_FOUND
    ERROR_VARIABLE cython_version
    ERROR_STRIP_TRAILING_WHITESPACE
    OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT CYTHON_FOUND EQUAL 0)
    message(FATAL_ERROR "Could not find cython")
else()
    set(CYTHON_FOUND True)
endif()
string(REGEX MATCH "[0-9].[0-9]*.[0-9]*" CYTHON_VERSION "${cython_version}")
message(STATUS "Found cython ${CYTHON_VERSION}")

# look for numpy
execute_process(
    COMMAND
    ${Python3_EXECUTABLE} -c "import numpy; print(numpy.__version__)"
    OUTPUT_VARIABLE NUMPY_VERSION
    RESULT_VARIABLE NUMPY_FOUND
    ERROR_VARIABLE numpy_error
    ERROR_STRIP_TRAILING_WHITESPACE
    OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT NUMPY_FOUND EQUAL 0)
    message(FATAL_ERROR "Could not find numpy")
else()
    set(NUMPY_FOUND True)
endif()
message(STATUS "Found numpy ${NUMPY_VERSION}")

# get numpy include directory
execute_process(
    COMMAND
    ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIRS
    RESULT_VARIABLE numpy_include_ok
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)
if(NOT numpy_include_ok EQUAL 0 OR NOT NUMPY_INCLUDE_DIRS)
    message(FATAL_ERROR "Could not find numpy includes")
endif()
message(STATUS "Found numpy includes at ${NUMPY_INCLUDE_DIRS}")

# add command to cythonize results
add_custom_command(
    OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/pyssht.c"
    COMMAND ${Python3_EXECUTABLE} -m cython
    -I "${CMAKE_CURRENT_SOURCE_DIR}/include/c/"
    -3
    -o "${CMAKE_CURRENT_BINARY_DIR}/pyssht.c"
    ${CMAKE_CURRENT_SOURCE_DIR}/src/python/pyssht.pyx
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/include/c/ssht.h"
    MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/src/python/pyssht.pyx
    COMMENT "Cythonizing src/python/pyssht.pyx"
    )
Python3_add_library(pyssht MODULE "${CMAKE_CURRENT_BINARY_DIR}/pyssht.c")
set_target_properties(pyssht PROPERTIES LIBRARY_OUTPUT_NAME ssht)
target_link_libraries(pyssht PUBLIC ssht)
target_include_directories(pyssht
    PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/include/c" ${NUMPY_INCLUDE_DIRS})
