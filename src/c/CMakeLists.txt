add_library(ssht STATIC ssht_adjoint.c ssht_core.c ssht_dl.c ssht_sampling.c)
target_link_libraries(ssht PUBLIC ${FFTW3_TARGET} m)
target_include_directories(
  ssht
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
         $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}>
         $<INSTALL_INTERFACE:include>)
set_target_properties(ssht PROPERTIES C_STANDARD 99)
configure_file(ssht_version.in.h ${CMAKE_CURRENT_BINARY_DIR}/ssht_version.h)

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
    FILES ssht_adjoint.h
          ssht_dl.h
          ssht_sampling.h
          ssht_core.h
          ssht_error.h
          ssht_types.h
          ${CMAKE_BINARY_DIR}/ssht_version.h
          ssht.h
    DESTINATION include/ssht)
endif()
