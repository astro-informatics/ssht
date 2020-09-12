# SSHT: Fast spin spherical harmonic transforms

[docs-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-url]: astro-informatics.github.io/ssht/

[![][docs-img]][docs-url]
![Bintray](https://img.shields.io/bintray/v/mdavezac/AstroFizz/ssht:AstroFizz?label=bintray%20-%20C%20package)
[![PyPI version](https://badge.fury.io/py/pyssht.svg)](https://badge.fury.io/py/pyssht)

## DESCRIPTION

The SSHT code provides functionality to perform fast and exact spin spherical
harmonic transforms.

## VERSION
Release 1.3.1, Sept 20

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

The C package can be installed with [CMake](https://cmake.org):

```cmake
git clone http://astro-informatics.github.io/ssht/ -b main
mkdir ssht/build && cd ssht/build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local ..
make
make install
```

The above will also download [FFTW](http://www.fftw.org/), if necessary.

Instructions for installing the fortran package can be found in docs/index.html.
