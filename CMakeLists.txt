cmake_minimum_required(VERSION 3.12)

project(Ssht
    VERSION "1.1.2"
    DESCRIPTION "Fast and exact spin spherical harmonic transforms"
    HOMEPAGE_URL "http://astro-informatics.github.io/ssht/"
    LANGUAGES C CXX)


set(Ssht_BUILD "CMake")
set(CMAKE_C_STANDARD 99)

option(openmp         "Enable OpenMP"                  ON)
option(fPIC           "Enable fPIC (only for gcc)"     ON)
option(tests          "Enable testing"                 ON)
option(python         "Creates python package"         OFF)
option(conan_fftw     "Download fftw using conan"      OFF)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

if(conan_fftw)
    if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
        message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
        file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake"
            "${CMAKE_BINARY_DIR}/conan.cmake")
    endif()
    include(${CMAKE_BINARY_DIR}/conan.cmake)

    if(NOT CONAN_OPTIONS)
        set(CONAN_OPTIONS
            "fftw:fPIC=True" "fftw:shared=False" "fftw:precision=double"
            "fftw:openmp=False")
    endif()
    if(NOT CONAN_BUILD)
        set(CONAN_BUILD "missing")
    endif()
    set(_GLIBCXX_USE_CXX11_ABI TRUE)

    conan_check(REQUIRED)
    conan_cmake_run(REQUIRES fftw/3.3.8
        BASIC_SETUP
        OPTIONS "${CONAN_OPTIONS}"
        KEEP_RPATHS
        CMAKE_TARGETS
        NO_OUTPUT_DIRS
        BUILD ${CONAN_BUILD}
        # import headers to help vim find them
        IMPORTS "include, * -> conan/include")
endif()

find_package(FFTW3 REQUIRED COMPONENT SERIAL DOUBLE)
if(openmp)
    find_package(OpenMP QUIET)
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
        # get numpy include directory
        execute_process(
            COMMAND
            ${Python3_EXECUTABLE} -c "import skbuild; print(skbuild.__file__)"
            OUTPUT_VARIABLE SKBUILD_LOCATION
            RESULT_VARIABLE SKBUILD_FOUND
            ERROR_QUIET
            OUTPUT_STRIP_TRAILING_WHITESPACE)
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
