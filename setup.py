from skbuild import setup

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
    cmake_args=[
        "-Dpython:BOOL=ON",
        "-Dopenmp:BOOL=OFF",
        "-Dtests:BOOL=OFF",
        "-DCMAKE_OSX_DEPLOYMENT_TARGET:STRING=10.9",
    ],
    cmake_languages=("C",),
    license="GPL-2",
    packages=["pyssht"],
)
