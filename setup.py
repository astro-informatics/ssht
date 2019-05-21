from sys import platform

from skbuild import setup

cmake_args = ["-Dpython:BOOL=ON", "-Dopenmp:BOOL=OFF", "-Dtests:BOOL=OFF"]
if platform.lower() == "darwin":
    cmake_args.append("-DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.9")

setup(
    name="pyssht",
    version="2.0",
    author="Jason McEwen",
    install_requires=["numpy", "cython"],
    extras_require={
        "dev": ["setuptools", "wheel", "scikit-build", "cmake", "ninja", "cython"]
    },
    description="Fast spin spherical transforms",
    url="http://astro-informatics.github.io/ssht/",
    package_dir={"pyssht": "src/pyssht"},
    package_data={"pyssht": ["SSHT_Python_Documentation.md"]},
    cmake_args=cmake_args,
    cmake_languages=("C",),
    license="GPL-2",
    packages=["pyssht"],
)
