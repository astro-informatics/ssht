if(NOT EXISTS "${CMAKE_BINARY_DIR}/conan.cmake")
  message(
    STATUS
      "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
  file(DOWNLOAD "https://github.com/conan-io/cmake-conan/raw/v0.16.1/conan.cmake"
       "${CMAKE_BINARY_DIR}/conan.cmake" TLS_VERIFY ON)
endif()
include(${CMAKE_BINARY_DIR}/conan.cmake)
conan_cmake_run(CONANFILE ${PROJECT_SOURCE_DIR}/conanfile.txt
                BASIC_SETUP CMAKE_TARGETS
                BUILD missing)
conan_cmake_autodetect(settings)
conan_cmake_install(PATH_OR_REFERENCE ${PROJECT_SOURCE_DIR}/conanfile.txt
                    BUILD missing
                    REMOTE conan-center
                    SETTINGS ${settings})