#ifndef SSHT_TEST_UTILITIES
#define SSHT_TEST_UTILITIES
#include <complex.h>

double ran2_dp(int idum);
void gen_flm_complex(complex double *flm, int L, int spin, int seed);
void gen_flm_real(complex double *flm, int L, int seed);
void gen_lb_flm_complex(complex double *flm, int L_zero, int L, int spin,
                        int seed);
void gen_lb_flm_real(complex double *flm, int L0, int L, int seed);
int max(int a, int b);
#endif
