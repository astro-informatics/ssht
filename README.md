# SSHT: Fast spin spherical harmonic transforms

[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: https://astro-informatics.github.io/ssht/
[conan-img]: https://img.shields.io/badge/ConanCenter-C%20Package-red.svg
[conan-url]: https://conan.io/center/ssht
[pypi-img]: https://badge.fury.io/py/pyssht.svg
[pypi-url]: https://badge.fury.io/py/pyssht
[codefactor-img]: https://www.codefactor.io/repository/github/astro-informatics/ssht/badge/main
[codefactor-url]: https://www.codefactor.io/repository/github/astro-informatics/ssht/overview/main

[![][docs-img]][docs-url]
[![][conan-img]][conan-url]
[![][pypi-img]][pypi-url]
[![][codefactor-img]][codefactor-url]
![CMake Build](https://github.com/astro-informatics/ssht/workflows/CMake%20Build/badge.svg)
![Python Build](https://github.com/astro-informatics/ssht/workflows/Python%20Build/badge.svg)

## DESCRIPTION

The **SSHT** code provides functionality to perform fast and exact
spin spherical harmonic transforms based on the sampling theorem on the
sphere derived in <a href="http://www.jasonmcewen.org/publication/mcewen-fssht/">McEwen & Wiaux (2011)</a>.

**SSHT** can also interface with [ducc0](https://pypi.org/project/ducc0/) and
use it as a backend for spherical harmonic transforms and rotations.


## INSTALLATION

 The python package, **pyssht**, is available on <a href="https://pypi.org/project/pyssht/">pypi</a> and can be installed with: 
 
 ```bash
 pip install pyssht
 ```

The C package can be installed with [CMake](https://cmake.org) and
[conan](https://docs.conan.io/en/latest/howtos/other_languages_package_manager/python.html).

Both can be installed using pip:

```bash
pip install conan cmake
```

Then **SSHT** can be compiled with:

```bash
git clone https://github.com/astro-informatics/ssht.git
mkdir ssht/build && cd ssht/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -Dconan_deps=True ..
make
make install
```

The above will also download [FFTW](http://www.fftw.org/), if necessary.


Further instructions for installing **SSHT** Matlab and Fortran versions are available [here](https://astro-informatics.github.io/ssht/).


## DOCUMENTATION

Further documentation is available [here](https://astro-informatics.github.io/ssht/).

Usage for the python package is also given in the package docstring.


## REFERENCING

If you use **SSHT** for work that results in publication, please reference <a href="https://github.com/astro-informatics/ssht">https://github.com/astro-informatics/ssht/</a> and cite our related academic papers:

- J. D. McEwen and Y. Wiaux, [A novel sampling theorem on the sphere](http://www.jasonmcewen.org/publication/mcewen-fssht/), IEEE Trans. Sig. Proc., 59(12):5876-5887, 2011 [(arXiv:1110.6298)](https://arxiv.org/abs/1110.6298).
- J. D. McEwen, G. Puy, J.-Ph. Thiran, P. Vandergheynst, D. Van De Ville, and Y. Wiaux, [Sparse image reconstruction on the sphere: implications of a new sampling theorem](http://www.jasonmcewen.org/publication/mcewen-css-2/). IEEE Trans. Image Proc., 22(6):2275-2285, 2013 [(arXiv:1205.1013)](https://arxiv.org/abs/1205.1013).


## LICENSE

SSHT is released under the GPL-3 license.  For further details see 
[LICENSE](https://github.com/astro-informatics/ssht/blob/main/LICENSE).


## AUTHORS

SSHT was initially written by [Jason McEwen](http://www.jasonmcewen.org/) but significant contributors have since been made by a number of <a href="https://github.com/astro-informatics/ssht/graphs/contributors">others</a>.
