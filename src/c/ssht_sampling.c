/*! 
 * \file ssht_sampling.c
 * Functionality to define sample positions for various algorithms,
 * to compute weights and to convert 1D and 2D harmonic indices.
 *
 * \author Jason McEwen
 */


#include <complex.h>
#include <math.h>
#include "ssht_types.h"


//============================================================================
// Sampling weights
//============================================================================


/*!
 * Compute weights for toroidal extension.
 *
 * \param[in]  p Integer index to compute weight for.
 * \param[out] w Corresponding weight.
 * \retval none
 *
 * \author Jason McEwen
 */
complex double ssht_sampling_weight_mw(int p) {

  if (p == 1) {
    return I * SSHT_PION2;
  }
  else if (p == -1) {
     return - I * SSHT_PION2;
  }
  else if (p % 2 == 0) {
    return 2.0 / (1.0 - p*p);
  }
  else {
    return 0.0;
  }
  
}
  

//============================================================================
// Sampling relations
//============================================================================


/*!
 * Convert theta index to angle for McEwen and Wiaux sampling.
 *
 * \note
 *  - t ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
 *
 * \param[in] t Theta index.
 * \param[in] L Harmonic band-limit.
 * \retval theta Theta angle.
 *
 * \author Jason McEwen
 */
double ssht_sampling_mw_t2theta(int t, int L) {

  return (2.0*t + 1.0) * SSHT_PI / (2.0*L - 1.0);

}


/*!
 * Convert phi index to angle for McEwen and Wiaux sampling.
 *
 * \note
 *  - p ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
 *
 * \param[in] p Phi index.
 * \param[in] L Harmonic band-limit.
 * \retval phi Phi angle.
 *
 * \author Jason McEwen
 */
double ssht_sampling_mw_p2phi(int p, int L) {

  return 2.0 * p * SSHT_PI / (2.0*L - 1.0);

}


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
 * \author Jason McEwen
 */
inline void ssht_sampling_elm2ind(int *ind, int el, int m) {

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
 * \author Jason McEwen
 */
inline void ssht_sampling_ind2elm(int *el, int *m, int ind) {

  *el = sqrt(ind);
  *m = ind - (*el)*(*el) - (*el);

}

     
