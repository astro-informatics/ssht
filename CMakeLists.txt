cmake_minimum_required(VERSION 3.7)

project(Ssht
    VERSION "1.1.2"
    DESCRIPTION "Fast and exact spin spherical harmonic transforms"
    HOMEPAGE_URL "http://astro-informatics.github.io/ssht/"
    LANGUAGES C)


set(Ssht_BUILD "CMake")
set(CMAKE_C_STANDARD 99)

option(openmp         "Enable OpenMP"                  ON)
option(fPIC           "Enable fPIC (only for gcc)"     ON)
option(tests          "Enable testing"                 ON)
option(python         "Creates python package"         OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(FFTW3 REQUIRED COMPONENT SERIAL DOUBLE)
find_package(OpenMP QUIET)

add_subdirectory(src/c)

if(python)
    add_subdirectory(src/python)
endif()
