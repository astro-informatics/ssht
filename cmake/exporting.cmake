# Exports ssht so other packages can access it
export(
    TARGETS ssht
    FILE "${PROJECT_BINARY_DIR}/sshtTargets.cmake"
    NAMESPACE ssht::)

# Avoids creating an entry in the cmake registry.
if(NOT NOEXPORT)
    export(PACKAGE ssht)
endif()

set(INCLUDE_INSTALL_DIR include/)
include(CMakePackageConfigHelpers)
configure_package_config_file(
    cmake/sshtConfig.in.cmake "${PROJECT_BINARY_DIR}/sshtConfig.cmake"
    INSTALL_DESTINATION lib/cmake/ssht
    PATH_VARS INCLUDE_INSTALL_DIR)
write_basic_package_version_file(
    sshtConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion)

if(NOT CONAN_EXPORTED)
    install(FILES "${PROJECT_BINARY_DIR}/sshtConfig.cmake"
                  "${PROJECT_BINARY_DIR}/sshtConfigVersion.cmake"
            DESTINATION lib/cmake/ssht)
endif()

install(
    EXPORT sshtTargets
    NAMESPACE ssht::
    DESTINATION lib/cmake/ssht)
