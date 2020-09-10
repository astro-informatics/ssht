cmake_minimum_required(VERSION 3.12)

project(
  Ssht
  VERSION "1.3.1"
  DESCRIPTION "Fast and exact spin spherical harmonic transforms"
  HOMEPAGE_URL "http://astro-informatics.github.io/ssht/"
  LANGUAGES C)

set(Ssht_BUILD "CMake")
set(CMAKE_C_STANDARD 99)

option(tests "Enable testing" ON)
option(python "Creates python package" OFF)
option(fPIC "Compile with fPIC" ON)
option(conan_deps "Download dependencies using conan" OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()
if(conan_deps)
    include("conan_dependencies")
else()
    find_package(FFTW3 REQUIRED COMPONENTS SERIAL DOUBLE)
endif()

add_subdirectory(src/c)
if(NOT python)
  include("exporting")
else()
  set(OLD_CMAKE_FIND_FRAMEWORK "${CMAKE_FIND_FRAMEWORK}")
  set(CMAKE_FIND_FRAMEWORK LAST)
  find_package(Python3 REQUIRED COMPONENTS Development Interpreter)
  set(CMAKE_FIND_FRAMEWORK "${OLD_CMAKE_FIND_FRAMEWORK}")

  if(NOT SKBUILD)
    execute_process(
      COMMAND ${Python3_EXECUTABLE} -c "import skbuild; print(skbuild.__file__)"
      OUTPUT_VARIABLE SKBUILD_LOCATION
      RESULT_VARIABLE SKBUILD_FOUND
      ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NOT SKBUILD_FOUND EQUAL 0 OR NOT SKBUILD_LOCATION)
      message(FATAL_ERROR "Could not find numpy includes")
    else()
      set(SKBUILD_FOUND True)
      get_filename_component(SKBUILD_LOCATION "${SKBUILD_LOCATION}" DIRECTORY)
    endif()
    message(STATUS "Found skbuild at ${SKBUILD_LOCATION}")
    list(APPEND CMAKE_MODULE_PATH "${SKBUILD_LOCATION}/resources/cmake")
  endif()

  find_package(PythonExtensions REQUIRED)
  find_package(Cython REQUIRED)
  find_package(NumPy REQUIRED)

  add_subdirectory(src/pyssht)
endif()
