from skbuild import setup

cmake_args = [
    "-Dpython:BOOL=ON",
    "-Dopenmp:BOOL=OFF",
    "-Dtests:BOOL=OFF",
    "-Dconan_fftw=ON",
]

build_requirements = [
    "setuptools",
    "wheel",
    "scikit-build",
    "cmake>=3.10",
    "ninja",
    "cython",
    "conan",
    "pip!=20.0.0,!=20.0.1",
]

setup(
    name="pyssht",
    version="2.0",
    author="Jason McEwen",
    install_requires=["numpy", "cython", "scipy"],
    extras_require={"build": build_requirements, "dev": build_requirements},
    description="Fast spin spherical transforms",
    url="http://astro-informatics.github.io/ssht/",
    package_dir={"pyssht": "src/pyssht"},
    package_data={"pyssht": ["SSHT_Python_Documentation.md"]},
    cmake_args=cmake_args,
    cmake_languages=("C",),
    license="GPL-2",
    packages=["pyssht"],
)
