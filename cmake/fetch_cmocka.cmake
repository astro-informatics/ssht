include(FetchContent)

FetchContent_Declare(
  CMocka
  GIT_REPOSITORY https://gitlab.com/cmocka/cmocka.git
  GIT_TAG cmocka-1.1.5)

if(NOT cmocka_POPULATED)
  # Fetch the content using previously declared details
  FetchContent_Populate("cmocka")
  file(READ "${cmocka_SOURCE_DIR}/CMakeLists.txt" patch_cmocka )
  string(REPLACE "add_subdirectory(doc)" "" patch_cmocka ${patch_cmocka})
  file(WRITE "${cmocka_SOURCE_DIR}/CMakeLists.txt" ${patch_cmocka} )

  add_subdirectory(${cmocka_SOURCE_DIR} ${cmocka_BINARY_DIR})
endif()

find_file(CSTDLIB_HEADERS NAMES setjmp.h PATHS /usr/include)