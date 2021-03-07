from conans import CMake, ConanFile


class SshtConan(ConanFile):
    name = "ssht"
    version = "1.3.5"
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
        "include/ssht/*.h",
        "tests/CMakeLists.txt",
        "tests/*.c",
        "tests/*.h",
        "CMakeLists.txt",
        "cmake/*.cmake",
    ]

    def configure(self):
        if self.settings.compiler == "Visual Studio":
            del self.options.fPIC
        self.options["fftw"].fPIC = self.options.fPIC
        del self.settings.compiler.libcxx

    @property
    def cmake(self):
        if not hasattr(self, "_cmake"):
            self._cmake = CMake(self)
            self._cmake.definitions["tests"] = True
            self._cmake.definitions["conan_deps"] = True
            self._cmake.definitions["python"] = False
            self._cmake.definitions["fPIC"] = self.options.fPIC
            self._cmake.configure(build_folder="build")
        return self._cmake

    def build(self):
        from pathlib import Path

        path = Path(self.source_folder)
        build = Path(self.source_folder) / "build"
        build.mkdir(exist_ok=True)
        (path / "conanbuildinfo.cmake").rename(path / "build" / "conanbuildinfo.cmake")
        self.cmake.build()
        self.cmake.test()

    def package(self):
        self.cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["ssht"]
