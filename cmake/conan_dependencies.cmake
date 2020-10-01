if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
  message(
    STATUS
      "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
  file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.15/conan.cmake"
       "${CMAKE_BINARY_DIR}/conan.cmake" TLS_VERIFY ON)
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)

if(fPIC)
  set(fpic_value "True")
else()
  set(fpic_value "False")
endif()
if(NOT CONAN_OPTIONS)
  set(CONAN_OPTIONS "fftw:shared=False" "fftw:precision=double"
                    "fftw:openmp=False")
  if(NOT WIN32)
    list(APPEND CONAN_OPTIONS "fftw:fPIC=${fpic_value}")
  endif()
endif()
if(NOT CONAN_BUILD)
    set(CONAN_BUILD "missing")
endif()

conan_check(REQUIRED)
conan_cmake_run(
  REQUIRES
  "fftw/3.3.8"
  BASIC_SETUP
  OPTIONS
  "${CONAN_OPTIONS}"
  KEEP_RPATHS
  CMAKE_TARGETS
  NO_OUTPUT_DIRS
  BUILD
  ${CONAN_BUILD})
