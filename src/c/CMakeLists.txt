add_library(ssht STATIC ssht_adjoint.c ssht_core.c ssht_dl.c ssht_sampling.c)
target_link_libraries(ssht PUBLIC fftw3 m)
target_include_directories(ssht PUBLIC $<INSTALL_INTERFACE:include>)

if(fPIC)
    set_target_properties(ssht PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()

if(NOT python)
    install(
        TARGETS ssht EXPORT SshtTargets ARCHIVE DESTINATION lib PUBLIC_HEADER)
    install(FILES
        ssht_adjoint.h ssht_dl.h ssht_sampling.h ssht_core.h ssht_error.h
        ssht_types.h ssht.h
        DESTINATION include/ssht)

    add_executable(ssht_about ssht_about.c)
    target_compile_definitions(ssht_about
        PRIVATE SSHT_VERSION="${Ssht_VERSION}" SSHT_BUILD="${Ssht_BUILD}")
    set_target_properties(ssht_about
        PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()

if(tests AND NOT python)
    add_executable(ssht_test ssht_test.c)
    target_link_libraries(ssht_test ssht)
    set_target_properties(ssht_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
    add_test(NAME so3_test COMMAND so3_test)
endif()
