add_cython_target(cython_sources pyssht.pyx C PY3)
python3_add_library(pyssht MODULE ${cython_sources})
target_link_libraries(pyssht PUBLIC ssht)
target_include_directories(pyssht PUBLIC "${PROJECT_SOURCE_DIR}/src/c"
                                         ${NumPy_INCLUDE_DIRS})
if(WIN32)
    target_compile_definitions(pyssht PUBLIC MS_WIN64)
endif()
install(TARGETS pyssht LIBRARY DESTINATION src/pyssht)
