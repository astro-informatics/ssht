add_library(ssht STATIC ssht_adjoint.c ssht_core.c ssht_dl.c ssht_sampling.c)
target_link_libraries(ssht PUBLIC FFTW3::FFTW3 ${MATH_LIBRARY})
target_include_directories(
    ssht
    PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
           $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
           $<INSTALL_INTERFACE:include>)
set_target_properties(ssht PROPERTIES ARCHIVE_OUTPUT_DIRECTORY
                                      ${PROJECT_BINARY_DIR}/lib)
target_compile_features(ssht PUBLIC c_std_99)

configure_file(${PROJECT_SOURCE_DIR}/include/ssht/ssht_version.in.h
               ${PROJECT_BINARY_DIR}/include/ssht/ssht_version.h)

if(NOT SKBUILD)
    install(
        TARGETS ssht
        EXPORT sshtTargets
        ARCHIVE DESTINATION lib
        PUBLIC_HEADER)
    install(
        FILES ${PROJECT_SOURCE_DIR}/include/ssht/ssht_adjoint.h
              ${PROJECT_SOURCE_DIR}/include/ssht/ssht_dl.h
              ${PROJECT_SOURCE_DIR}/include/ssht/ssht_sampling.h
              ${PROJECT_SOURCE_DIR}/include/ssht/ssht_core.h
              ${PROJECT_SOURCE_DIR}/include/ssht/ssht_error.h
              ${PROJECT_SOURCE_DIR}/include/ssht/ssht_types.h
              ${PROJECT_BINARY_DIR}/include/ssht/ssht_version.h
              ${PROJECT_SOURCE_DIR}/include/ssht/ssht.h
        DESTINATION include/ssht)
endif()
