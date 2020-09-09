cmake_minimum_required(VERSION 2.8.12)
project(PackageTest C)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

find_package(Ssht)
add_executable(example example.c)
target_link_libraries(example ${CONAN_LIBS} ssht)

add_test(NAME example
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/bin
         COMMAND example)
