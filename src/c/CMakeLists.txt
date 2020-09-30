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

  add_executable(ssht_about ssht_about.c)
  target_include_directories(ssht_about PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
  set_source_files_properties(
    ssht_about.c PROPERTIES OBJECT_DEPENDS
                            ${CMAKE_CURRENT_BINARY_DIR}/ssht_version.h)
  set_target_properties(ssht_about PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                              ${PROJECT_BINARY_DIR})
  target_compile_definitions(ssht_about PRIVATE BUILT_WITH_CMAKE)
endif()

if(tests AND NOT python)
  add_executable(ssht_test ssht_test.c)
  target_link_libraries(ssht_test ssht)
  set_target_properties(ssht_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                             ${PROJECT_BINARY_DIR})
  add_test(NAME ssht_test COMMAND ssht_test)
endif()
