if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.15/conan.cmake"
                 "${CMAKE_BINARY_DIR}/conan.cmake" 
                 TLS_VERIFY ON)
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)

if(NOT CONAN_OPTIONS)
    set(CONAN_OPTIONS "fftw:shared=False" "fftw:precision=double" "fftw:openmp=False")
    if(fPIC)
        list(APPEND CONAN_OPTIONS "fftw:fPIC=True")
    else()
        list(APPEND CONAN_OPTIONS "fftw:fPIC=False")
    endif()
endif()
if(NOT CONAN_DEPS)
    set(CONAN_DEPS "fftw/3.3.8")
endif()
if(NOT CONAN_BUILD)
    set(CONAN_BUILD "missing")
endif()

conan_check(REQUIRED)
conan_cmake_run(REQUIRES ${CONAN_DEPS}
    BASIC_SETUP
    OPTIONS "${CONAN_OPTIONS}"
    KEEP_RPATHS
    CMAKE_TARGETS
    NO_OUTPUT_DIRS
    BUILD ${CONAN_BUILD})
add_library(fftw3 ALIAS CONAN_PKG::fftw)
