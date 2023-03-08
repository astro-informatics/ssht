from conan import ConanFile
from conan.tools.cmake import CMake, CMakeDeps, CMakeToolchain, cmake_layout


class Ssht(ConanFile):
    name = "ssht"
    version = "1.5.1"
    description = "Fast and exact spin spherical harmonic transforms."
    license = "GPL-3"
    author = "Jason McEwen"
    homepage = "https://astro-informatics.github.io/ssht/"
    url = "https://github.com/astro-informatics/ssht"
    settings = "os", "arch", "compiler", "build_type"
    requires = "fftw/[>=3.3]"
    tool_requires = "cmake/[>=3.18]"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = {
        "shared": False,
        "fPIC": True,
        "fftw*:shared": False,
        "fftw*:precision": "double",
    }

    def generate(self):
        tc = CMakeToolchain(self)
        tc.generate()
        deps = CMakeDeps(self)
        deps.generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def layout(self):
        cmake_layout(self)
