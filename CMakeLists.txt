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
if(tests)
    enable_testing()
endif()
if(conan_deps)
    include("conan_dependencies")
    set(FFTW3_TARGET CONAN_PKG::fftw)
else()
    find_package(FFTW3 REQUIRED COMPONENTS SERIAL DOUBLE)
    set(FFTW3_TARGET fftw3)
endif()

add_subdirectory(src/c)

if(NOT python)
  include("exporting")
else()
  include("python-setup")
  find_python()
  setup_skbuild()
  find_package(PythonExtensions REQUIRED)
  find_package(Cython REQUIRED)
  find_package(NumPy REQUIRED)
  add_subdirectory(src/pyssht)
endif()
