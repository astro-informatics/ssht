add_library(ssht STATIC ssht_adjoint.c ssht_core.c ssht_dl.c ssht_sampling.c)
target_link_libraries(ssht ${FFTW3_DOUBLE_SERIAL_LIBRARY})
target_link_libraries(ssht m)
if(OpenMP_C_FOUND AND openmp)
    target_link_libraries(ssht OpenMP::OpenMP_C)
endif()

target_include_directories(ssht PUBLIC ${FFTW3_INCLUDE_DIRS})
if(fPIC)
    set_target_properties(ssht
        PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS} -fPIC")
else()
    set_target_properties(ssht PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS}")
endif()

if(NOT python)
    install(
        TARGETS ssht EXPORT SshtTargets ARCHIVE DESTINATION lib PUBLIC_HEADER)
    install(FILES
        ssht_adjoint.h ssht_dl.h ssht_sampling.h ssht_core.h ssht_error.h
        ssht_types.h ssht.h
        DESTINATION include/ssht)

    #building about and test executables
    add_executable(ssht_about ssht_about.c)
    target_compile_definitions(ssht_about
        PRIVATE SSHT_VERSION="${Ssht_VERSION}" SSHT_BUILD="${Ssht_BUILD}")
    set_target_properties(ssht_about
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

if(tests AND NOT python)
    add_executable(ssht_test ssht_test.c)
    target_include_directories(ssht_test PUBLIC ${FFTW3_INCLUDE_DIRS})
    target_link_libraries(ssht_test ${FFTW3_DOUBLE_SERIAL_LIBRARY})
    set_target_properties(ssht_test
        PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS}")
    if(OpenMP_C_FOUND AND openmp)
        target_link_libraries(ssht_test OpenMP::OpenMP_C)
    endif()
    target_link_libraries(ssht_test ssht)
    set_target_properties(ssht_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()
