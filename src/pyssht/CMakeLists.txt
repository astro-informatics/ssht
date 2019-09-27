# add command to cythonize results
add_cython_target(cython_sources pyssht.pyx C PY3)

Python3_add_library(pyssht MODULE ${cython_sources})
target_link_libraries(pyssht PUBLIC ssht)
target_include_directories(pyssht
    PUBLIC "${PROJECT_SOURCE_DIR}/include/c" ${NumPy_INCLUDE_DIRS})
install(TARGETS pyssht LIBRARY DESTINATION src/pyssht)
