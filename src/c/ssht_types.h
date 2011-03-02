/*! \mainpage SSHT: Spin spherical harmonic transforms
 *
 * The SSHT code provides functionality to perform fast and exact spin
 * spherical harmonic transforms based on the sampling theorem on the
 * sphere derived in our paper: <i>A novel sampling theorem on the
 * sphere</i> (<a href="http://arxiv.org/abs/XXX.XXX">ArXiv</a>|
 * <a href="http://dx.doi.org/10.1111/XXX">DOI</a>).
 *   
 * We document the C source code here.  For an example of usage, 
 * see the ssht_test.c program.
 * For installation instructions, see the general SSHT 
 * documentation available 
 * <a href="../../index_ssht.html">here</a>.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 * \version 0.1
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


#define SSHT_VERSION 0.1
#define SSHT_PROMPT "[ssht-0.1]"



#endif
