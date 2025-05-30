cmake_minimum_required(VERSION 3.24)
if(NOT SKBUILD_PROJECT_NAME)
    set(project_name Ssht)
else()
    set(project_name ${SKBUILD_PROJECT_NAME})
endif()
project(
    ${project_name}
    VERSION "1.5.3"
    DESCRIPTION "Fast and exact spin spherical harmonic transforms"
    HOMEPAGE_URL "http://astro-informatics.github.io/ssht/"
    LANGUAGES C)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()
include(CTest)

if(EXISTS "${PROJECT_BINARY_DIR}/FindFFTW3.cmake")
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_BINARY_DIR}")
else()
    list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake")
endif()
find_package(FFTW3 REQUIRED)
find_library(MATH_LIBRARY m)

add_subdirectory(src/c)
if(BUILD_TESTING)
    include("${PROJECT_SOURCE_DIR}/cmake/fetch_cmocka.cmake")
    add_subdirectory(tests)
endif()

if(SKBUILD)
    find_package(
        Python
        COMPONENTS Interpreter Development.Module NumPy
        REQUIRED)
    add_subdirectory(src/pyssht)
else()
    include("${PROJECT_SOURCE_DIR}/cmake/exporting.cmake")
endif()

# only run documentation if this is not a sub-project
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR AND NOT SKBUILD)
    find_package(Doxygen)
    if(DOXYGEN_FOUND)
        set(DOXYGEN_OPTIMIZE_OUTPUT_FOR_C "YES")
        set(DOXYGEN_EXTRACT_ALL "YES")
        set(DOXYGEN_FILE_PATTERNS *.h *.c)
        doxygen_add_docs(docs src/c include tests)
    endif()
endif()
