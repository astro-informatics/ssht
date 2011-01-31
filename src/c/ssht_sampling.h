
#ifndef SSHT_SAMPLING
#define SSHT_SAMPLING

#include <complex.h>

complex double ssht_sampling_weight_mw(int p);
double ssht_sampling_mw_t2theta(int t, int L);
extern inline void ssht_sampling_elm2ind(int *ind, int el, int m);
extern inline void ssht_sampling_ind2elm(int *el, int *m, int ind);

#endif
