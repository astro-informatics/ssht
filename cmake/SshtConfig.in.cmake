set(Ssht_VERSION "@PROJECT_VERSION@")

@PACKAGE_INIT@

if(NOT "@conan_deps@" AND NOT "@CONAN_EDITABLE_MODE@")
    include(CMakeFindDependencyMacro)
    find_dependency(FFTW3 REQUIRED COMPONENTS SERIAL DOUBLE)
endif()

set(Ssht_INCLUDE_DIRS "@PACKAGE_INCLUDE_INSTALL_DIR@")
set(Ssht_LIBRARIES ssht)

check_required_components(Ssht)
