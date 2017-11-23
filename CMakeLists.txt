cmake_minimum_required(VERSION 3.7)
option(openmp         "Enable OpenMP"                                   on)
option(fPIC           "Enable fPIC (only for gcc)"                      on)
option(tests          "Enable testing"                                  on)

project(ssht C)


set(Ssht_VERSION "1.1.2")
set(Ssht_BUILD "CMake")
set(CMAKE_C_STANDARD 99)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(FFTW3 REQUIRED COMPONENT SERIAL DOUBLE)
find_package(OpenMP QUIET)

add_library(ssht STATIC
  src/c/ssht_adjoint.c src/c/ssht_core.c src/c/ssht_dl.c src/c/ssht_sampling.c)
if(fPIC)
  SET_TARGET_PROPERTIES(ssht PROPERTIES COMPILE_FLAGS "-fPIC")
endif()
target_link_libraries(ssht ${FFTW3_DOUBLE_SERIAL_LIBRARY})
target_link_libraries(ssht m)
if(OpenMP_C_FOUND AND openmp)
  target_link_libraries(ssht OpenMP::OpenMP_C)
endif()

target_include_directories(ssht PUBLIC ${FFTW3_INCLUDE_DIRS})
set_target_properties(ssht PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS}")

install(TARGETS ssht EXPORT SshtTargets ARCHIVE DESTINATION lib PUBLIC_HEADER)
install(FILES
  src/c/ssht_adjoint.h  src/c/ssht_dl.h       src/c/ssht_sampling.h
  src/c/ssht_core.h     src/c/ssht_error.h    src/c/ssht_types.h
  src/c/ssht.h
  DESTINATION include/ssht)

include("exporting")
#building about and test executables
add_executable(ssht_about src/c/ssht_about.c)
target_compile_definitions(ssht_about PRIVATE SSHT_VERSION="${Ssht_VERSION}" SSHT_BUILD="${Ssht_BUILD}")
set_target_properties(ssht_about PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
if(tests)
  add_executable(ssht_test src/c/ssht_test.c)
  target_include_directories(ssht_test PUBLIC ${FFTW3_INCLUDE_DIRS})
  target_link_libraries(ssht_test ${FFTW3_DOUBLE_SERIAL_LIBRARY})
  set_target_properties(ssht_test PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS}")
  if(OpenMP_C_FOUND AND openmp)
    target_link_libraries(ssht_test OpenMP::OpenMP_C)
  endif()
  target_link_libraries(ssht_test ssht)
  set_target_properties(ssht_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR})
endif()
