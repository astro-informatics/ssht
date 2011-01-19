/*! \file ssht_sampling.h
 *  Functionality to define sample positions for various algorithms,
 *  to compute weights and to convert 1D and 2D harmonic indices.
 *
 * \author Jason McEwen
 */

#ifndef SSHT_SAMPLING
#define SSHT_SAMPLING

#include <complex.h>


complex double ssht_sampling_weight_mw(int p);
extern inline void ssht_sampling_elm2ind(int *ind, int el, int m);
extern inline void ssht_sampling_ind2elm(int *el, int *m, int ind);


#endif
