from pathlib import Path

from skbuild import setup

cmake_args = [
    "-DBUILD_TESTING:BOOL=OFF",
    "-Dconan_deps=ON",
    "-DCMAKE_POSITION_INDEPENDENT_CODE=ON",
]

dev_requirements = [
    "setuptools",
    "wheel",
    "scikit-build",
    "cmake>=3.12",
    "ninja",
    "cython",
    "conan",
    "pip!=20.0.0,!=20.0.1",
    "pytest",
    "ducc0>=0.18",
]

long_description = (
    Path(__file__).parent / "src" / "pyssht" / "SSHT_Python_Documentation.md"
).read_text()

setup(
    name="pyssht",
    version="1.5.1",
    author=[
        "J. D. McEwen",
        "C. R. G. Wallis",
        "M. Buttner",
        "B. Leistedt",
        "Y. Wiaux",
    ],
    install_requires=["numpy", "scipy"],
    extras_require={"dev": dev_requirements, "ducc0": ["ducc0>=0.18"]},
    description="Fast spin spherical transforms",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://astro-informatics.github.io/ssht/",
    package_dir={"pyssht": "src/pyssht"},
    package_data={"pyssht": ["SSHT_Python_Documentation.md"]},
    cmake_args=cmake_args,
    cmake_languages=("C",),
    license="GPL-3",
    packages=["pyssht"],
)
