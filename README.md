# SSHT: Fast spin spherical harmonic transforms

[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: https://astro-informatics.github.io/ssht/
[bintray-img]: https://img.shields.io/bintray/v/mdavezac/AstroFizz/ssht:AstroFizz?label=C%20package
[bintray-url]: https://bintray.com/mdavezac/AstroFizz/ssht:AstroFizz/1.3.2:stable/link
[pypi-img]: https://badge.fury.io/py/pyssht.svg
[pypi-url]: https://badge.fury.io/py/pyssht

[![][docs-img]][docs-url]
[![][bintray-img]][bintray-url]
[![][pypi-img]][pypy-url]

## DESCRIPTION

The SSHT code provides functionality to perform fast and exact spin spherical
harmonic transforms.

## AUTHORS
- [J. D. McEwen](https://www.jasonmcewen.org)
- C. G. R. Wallis
- M. Buttner
- B. Leistedt
- Y. Wiaux

## REFERENCES
- J. D. McEwen and Y. Wiaux, A novel sampling theorem on the sphere, IEEE Trans. Sig. Proc., 59(12):5876-5887, 2011 [(arXiv:1110.6298)](https://arxiv.org/abs/1110.6298).
- J. D. McEwen, G. Puy, J.-Ph. Thiran, P. Vandergheynst, D. Van De Ville, and Y. Wiaux, Sparse image reconstruction on the sphere: implications of a new sampling theorem. IEEE Trans. Image Proc., 22(6):2275-2285, 2013 [(arXiv:1205.1013)](https://arxiv.org/abs/1205.1013).

## DOCUMENTATION
See docs/index.html

Usage for the python package is given in the package docstring.

## INSTALLATION
The python package can be installed with ``pip install pyssht``.

The C package can be installed with [CMake](https://cmake.org) and
[conan](https://docs.conan.io/en/latest/howtos/other_languages_package_manager/python.html):

Both can be installed using pip:

```bash
pip install conan cmake
```

Then ssht can be compiled with:

```bash
git clone http://astro-informatics.github.io/ssht/ -b main
mkdir ssht/build && cd ssht/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local -Dconan_deps=True ..
make
make install
```

The above will also download [FFTW](http://www.fftw.org/), if necessary.

Instructions for installing the fortran package can be found in docs/index.html.
