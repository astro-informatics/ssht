add_library(ssht STATIC ssht_adjoint.c ssht_core.c ssht_dl.c ssht_sampling.c)
target_link_libraries(ssht PUBLIC ${FFTW3_TARGET} m)
target_include_directories(
  ssht
  PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
         $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}/include>
         $<INSTALL_INTERFACE:include>)
set_target_properties(ssht PROPERTIES C_STANDARD 99 ARCHIVE_OUTPUT_DIRECTORY
                                                    ${PROJECT_BINARY_DIR}/lib)
configure_file(${PROJECT_SOURCE_DIR}/include/ssht/ssht_version.in.h
               ${PROJECT_BINARY_DIR}/include/ssht/ssht_version.h)

if(fPIC)
  set_target_properties(ssht PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
endif()

if(NOT python)
  install(
    TARGETS ssht
    EXPORT SshtTargets
    ARCHIVE DESTINATION lib
    PUBLIC_HEADER)
  install(
    FILES ${PROJECT_SOURCE_DIR}/include/ssht/ssht_adjoint.h
          ${PROJECT_SOURCE_DIR}/include/ssht/ssht_dl.h
          ${PROJECT_SOURCE_DIR}/include/ssht/ssht_sampling.h
          ${PROJECT_SOURCE_DIR}/include/ssht/ssht_core.h
          ${PROJECT_SOURCE_DIR}/include/ssht/ssht_error.h
          ${PROJECT_SOURCE_DIR}/include/ssht/ssht_types.h
          ${CMAKE_CURRENT_BINARY_DIR}/include/ssht/ssht_version.h
          ${PROJECT_SOURCE_DIR}/include/ssht/ssht.h
    DESTINATION include/ssht)
endif()
