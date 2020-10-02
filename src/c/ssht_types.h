// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details

/*! \mainpage SSHT C documentation
 *
 * The SSHT code provides functionality to perform fast and exact spin
 * spherical harmonic transforms based on the sampling theorem on the
 * sphere derived in our paper: <i><a
 * href="http://www.mrao.cam.ac.uk/~jdm57/publications.html#mcewen:fssht">A
 * novel sampling theorem on the sphere</a></i> (<a
 * href="http://arxiv.org/abs/1110.6298">ArXiv</a> | <a
 * href="http://dx.doi.org/10.1109/TSP.2011.2166394">DOI</a>).
 * Functionality is also provided to perform fast and exact adjoint
 * transforms based on the fast algorithms derived in our paper: <i><a
 * href="http://www.mrao.cam.ac.uk/~jdm57/publications.html#mcewen:css2">Efficient
 * and compressive sampling on the sphere</a></i>.
 *   
 * We document the C source code here.  For an example of usage, 
 * see the ssht_test.c program.
 * For installation instructions, see the general SSHT 
 * documentation available 
 * <a href="../../index.html">here</a>.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */


/*! \file ssht_types.h
 *  Types used in SSHT package.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#ifndef SSHT_TYPES
#define SSHT_TYPES

#define SSHT_PI    3.141592653589793238462643383279502884197
#define SSHT_PION2 1.570796326794896619231321691639751442099

#define SSHT_SQRT2 1.41421356237309504880168872420969807856967



#define SSHT_PROMPT "[ssht] "

#ifdef __cplusplus
#include <complex>
typedef std::complex<double> ssht_complex_double;
extern "C" {
#else
#include <tgmath.h>
typedef double complex ssht_complex_double;
#endif


#ifdef __cplusplus
}
#endif


#endif
