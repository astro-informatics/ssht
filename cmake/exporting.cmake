# Exports ssht so other packages can access it
export(TARGETS ssht FILE "${PROJECT_BINARY_DIR}/SshtTargets.cmake")

# Avoids creating an entry in the cmake registry.
if(NOT NOEXPORT)
  export(PACKAGE Ssht)
endif()

# First in binary dir
set(ALL_INCLUDE_DIRS "${PROJECT_SOURCE_DIR}")
configure_File(cmake/SshtConfig.in.cmake
  "${PROJECT_BINARY_DIR}/SshtConfig.cmake" @ONLY
)
configure_File(cmake/SshtConfigVersion.in.cmake
  "${PROJECT_BINARY_DIR}/SshtConfigVersion.cmake" @ONLY
)

# Then for installation tree
file(RELATIVE_PATH REL_INCLUDE_DIR
    "${CMAKE_INSTALL_PREFIX}/share/cmake/ssht"
    "${CMAKE_INSTALL_PREFIX}/include"
)
set(ALL_INCLUDE_DIRS "\${Ssht_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(cmake/SshtConfig.in.cmake
  "${PROJECT_BINARY_DIR}/CMakeFiles/SshtConfig.cmake" @ONLY
)

# Finally install all files
install(FILES
  "${PROJECT_BINARY_DIR}/CMakeFiles/SshtConfig.cmake"
  "${PROJECT_BINARY_DIR}/SshtConfigVersion.cmake"
    DESTINATION share/cmake/ssht
    COMPONENT dev
)

install(EXPORT SshtTargets DESTINATION share/cmake/ssht COMPONENT dev)
