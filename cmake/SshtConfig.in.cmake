get_filename_component(Ssht_CMAKE_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
message(STATUS "Linking to ssht package in ${Ssht_CMAKE_DIR}")
if(NOT TARGET ssht AND EXISTS "${Ssht_CMAKE_DIR}/SshtTargets.cmake")
  include("${Ssht_CMAKE_DIR}/SshtTargets.cmake")
endif()

set(Ssht_INCLUDE_DIRS "@ALL_INCLUDE_DIRS@")
set(Ssht_LIBRARIES ssht)
