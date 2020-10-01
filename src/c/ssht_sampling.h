// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SSHT_SAMPLING
#define SSHT_SAMPLING

#include <math.h>

#include "ssht_types.h"


#ifdef __cplusplus
extern "C"{
#endif

ssht_complex_double ssht_sampling_weight_mw(int p);
double ssht_sampling_weight_dh(double theta_t, int L);
void ssht_sampling_gl_thetas_weights(double *thetas, double *weights, int L);

double ssht_sampling_mw_t2theta(int t, int L);
double ssht_sampling_mw_p2phi(int p, int L);
int ssht_sampling_mw_n(int L);
int ssht_sampling_mw_ntheta(int L);
int ssht_sampling_mw_nphi(int L);

double ssht_sampling_mw_ss_t2theta(int t, int L);
double ssht_sampling_mw_ss_p2phi(int p, int L);
int ssht_sampling_mw_ss_n(int L);
int ssht_sampling_mw_ss_ntheta(int L);
int ssht_sampling_mw_ss_nphi(int L);

double ssht_sampling_dh_t2theta(int t, int L);
double ssht_sampling_dh_p2phi(int p, int L);
int ssht_sampling_dh_n(int L);
int ssht_sampling_dh_ntheta(int L);
int ssht_sampling_dh_nphi(int L);

double ssht_sampling_gl_p2phi(int p, int L);
int ssht_sampling_gl_n(int L);
int ssht_sampling_gl_ntheta(int L);
int ssht_sampling_gl_nphi(int L);


//============================================================================
// Harmonic index relations
//============================================================================


/*!
 * Convert (el,m) harmonic indices to 1D index used to access flm
 * array.
 *
 * \note Index ranges are as follows:  
 *  - el ranges from [0 .. L-1].
 *  - m ranges from [-el .. el].
 *  - ind ranges from [0 .. L**2-1].
 *
 * \param[out] ind 1D index to access flm array [output].
 * \param[in]  el  Harmonic index [input].
 * \param[in]  m   Azimuthal harmonic index [input].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
static inline void ssht_sampling_elm2ind(int *ind, int el, int m) {

  *ind = el * el + el + m;

}


/*!
 * Convert 1D index used to access flm array to (el,m) harmonic
 * indices.
 *
 * \note Index ranges are as follows:  
 *  - el ranges from [0 .. L-1].
 *  - m ranges from [-el .. el].
 *  - ind ranges from [0 .. L**2-1].
 *
 * \param[in]  ind 1D index to access flm array [output].
 * \param[out] el  Harmonic index [input].
 * \param[out] m   Azimuthal harmonic index [input].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
static inline void ssht_sampling_ind2elm(int *el, int *m, int ind) {

  *el = lrint(floor(sqrt(ind)));
  *m = ind - (*el)*(*el) - (*el);

}

#ifdef __cplusplus
}
#endif

#endif
