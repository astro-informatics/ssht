set(SSHT_VERSION "@PROJECT_VERSION@")

@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(FFTW3 REQUIRED)

set(SSHT_INCLUDE_DIRS "@PACKAGE_INCLUDE_INSTALL_DIR@")
set(SSHT_LIBRARIES ssht::ssht)

check_required_components(SSHT)
