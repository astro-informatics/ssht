set(Ssht_VERSION "@PROJECT_VERSION@")

@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(FFTW3 REQUIRED)

set(Ssht_INCLUDE_DIRS "@PACKAGE_INCLUDE_INSTALL_DIR@")
set(Ssht_LIBRARIES ssht)

check_required_components(Ssht)
