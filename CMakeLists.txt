cmake_minimum_required(VERSION 3.7)
project(ssht C)

set(Ssht_VERSION "1.1.2")
set(CMAKE_C_STANDARD 99)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(FFTW3 REQUIRED COMPONENT SERIAL DOUBLE)
find_package(OpenMP QUIET)

add_library(ssht STATIC
  src/c/ssht_adjoint.c src/c/ssht_core.c src/c/ssht_dl.c src/c/ssht_sampling.c)
target_link_libraries(ssht ${FFTW3_DOUBLE_SERIAL_LIBRARY})
if(OpenMP_C_FOUND)
  target_link_libraries(ssht OpenMP::OpenMP_C)
endif()

target_include_directories(ssht PUBLIC ${FFTW3_INCLUDE_DIRS})
set_target_properties(ssht PROPERTIES COMPILE_FLAGS "${FFTW3_DEFINITIONS}")

install(TARGETS ssht EXPORT SshtTargets ARCHIVE DESTINATION lib PUBLIC_HEADER)
install(FILES
  src/c/ssht_adjoint.h  src/c/ssht_dl.h       src/c/ssht_sampling.h
  src/c/ssht_core.h     src/c/ssht_error.h    src/c/ssht_types.h
  DESTINATION include/ssht)

include("exporting")
