# Exports ssht so other packages can access it
export(TARGETS ssht FILE "${PROJECT_BINARY_DIR}/SshtTargets.cmake")

# Avoids creating an entry in the cmake registry.
if(NOT NOEXPORT)
  export(PACKAGE Ssht)
endif()

set(INCLUDE_INSTALL_DIR include/)
include(CMakePackageConfigHelpers)
configure_package_config_file(cmake/SshtConfig.in.cmake
  "${PROJECT_BINARY_DIR}/SshtConfig.cmake"
  INSTALL_DESTINATION lib/cmake/Ssht
  PATH_VARS INCLUDE_INSTALL_DIR)
write_basic_package_version_file(
  SshtConfigVersion.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion)

if(NOT CONAN_EXPORTED)
  install(FILES
    "${PROJECT_BINARY_DIR}/SshtConfig.cmake"
    "${PROJECT_BINARY_DIR}/SshtConfigVersion.cmake"
    DESTINATION lib/cmake/ssht
    COMPONENT dev
  )
endif()

install(EXPORT SshtTargets DESTINATION lib/cmake/ssht COMPONENT dev)