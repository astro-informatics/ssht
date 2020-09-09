from conans import ConanFile, CMake


class SshtConan(ConanFile):
    name = "ssht"
    version = "1.3.0"
    license = "GPL-3.0"
    url = "https://github.com/astro-informatics/ssht"
    homepage = "https://github.com/astro-informatics/ssht"
    description = "Fast spin spherical harmonic transforms"
    settings = "os", "arch", "compiler", "build_type"
    topics = ("Physics", "Astrophysics", "Radio Interferometry")
    options = {"conan_fftw": [True, False]}
    default_options = "conan_fftw=True"
    generators = "cmake"
    exports_sources = [
        "src/c/*",
        "CMakeLists.txt",
        "src/pyssht/*",
        "cmake/*.cmake",
    ]

    def requirements(self):
        if self.options.conan_fftw:
            self.requires("fftw/3.3.8")

    def configured_cmake(self):
        cmake = CMake(self)
        cmake.definitions["CMAKE_POSITION_INDEPENDENT_CODE"] = True
        cmake.definitions["tests"] = False
        if self.options.conan_fftw:
            cmake.definitions["FFTW3_INCLUDE_DIR"] = ";".join(
                self.deps_cpp_info["fftw"].include_paths
            )
            cmake.definitions["FFTW3_LIBRARY_DIR"] = ";".join(
                self.deps_cpp_info["fftw"].lib_paths
            )
        cmake.configure(source_folder=".")
        return cmake

    def build(self):
        cmake = self.configured_cmake()
        cmake.build()

    def package(self):
        cmake = self.configured_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["ssht"]
