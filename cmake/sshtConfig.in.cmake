set(SSHT_VERSION "@PROJECT_VERSION@")

@PACKAGE_INIT@

include(CMakeFindDependencyMacro)
find_dependency(FFTW3 REQUIRED)

include("${CMAKE_CURRENT_LIST_DIR}/sshtTargets.cmake")
set(SSHT_LIBRARIES ssht::ssht)

check_required_components(SSHT)
