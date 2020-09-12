from conans import ConanFile, CMake


class SshtConan(ConanFile):
    name = "ssht"
    version = "1.3.1"
    license = "GPL-3.0"
    url = "https://github.com/astro-informatics/ssht"
    homepage = "https://github.com/astro-informatics/ssht"
    description = "Fast spin spherical harmonic transforms"
    settings = "os", "arch", "compiler", "build_type"
    topics = ("Physics", "Astrophysics", "Radio Interferometry")
    options = {"fPIC": [True, False]}
    default_options = {"fPIC": True}
    requires = "fftw/3.3.8"
    generators = "cmake"
    exports_sources = [
        "src/c/*",
        "CMakeLists.txt",
        "cmake/*.cmake",
    ]

    def configured_cmake(self):
        cmake = CMake(self)
        cmake.definitions["tests"] = True
        cmake.definitions["conan_deps"] = True
        cmake.definitions["python"] = False
        cmake.definitions["fPIC"] = self.options.fPIC
        cmake.configure(source_folder=".")
        return cmake

    def build(self):
        cmake = self.configured_cmake()
        cmake.build()
        cmake.test()

    def package(self):
        cmake = self.configured_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["ssht"]
