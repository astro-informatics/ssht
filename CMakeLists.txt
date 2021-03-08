cmake_minimum_required(VERSION 3.12)

project(
  Ssht
  VERSION "1.3.6rc1"
  DESCRIPTION "Fast and exact spin spherical harmonic transforms"
  HOMEPAGE_URL "http://astro-informatics.github.io/ssht/"
  LANGUAGES C)

option(tests "Enable testing" ON)
option(conan_deps "Download dependencies using conan" OFF)
if(NOT WIN32)
  option(fPIC "Compile with fPIC" ON)
endif()

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()
if(tests)
  enable_testing()
endif()
if(conan_deps OR CONAN_EDITABLE_MODE)
  include("conan_dependencies")
  set(FFTW3_TARGET CONAN_PKG::fftw)
else()
  find_package(FFTW3 REQUIRED COMPONENTS SERIAL DOUBLE)
  set(FFTW3_TARGET fftw3)
endif()

add_subdirectory(src/c)
if(tests)
  include("fetch_cmocka")
  add_subdirectory(tests)
endif()

if(NOT SKBUILD)
  include("exporting")
else()
  find_package(PythonExtensions REQUIRED)
  find_package(Cython REQUIRED)
  find_package(NumPy REQUIRED)
  add_subdirectory(src/pyssht)
endif()

# only run documentation if this is not a sub-project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  find_package(Doxygen)
  if(DOXYGEN_FOUND)
    set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C "YES")
    set(DOXYGEN_EXTRACT_ALL "YES")
    set(DOXYGEN_FILE_PATTERNS *.h *.c)
    doxygen_add_docs(docs src/c include tests)
  endif()
endif()
