
#ifndef SSHT_CORE
#define SSHT_CORE

#include <complex.h>

void ssht_core_mw_inverse_sov_sym(complex double *f, complex double *flm, 
				  int L, int spin, int verbosity);
void ssht_core_mw_inverse_sov_sym_real(double *f, complex double *flm, 
				       int L, int verbosity);
void ssht_core_mw_forward_sov_conv_sym(complex double *flm, complex double *f, 
				       int L, int spin, int verbosity);
void ssht_core_mw_forward_sov_conv_sym_real(complex double *flm, double *f, 
					    int L, int verbosity);


#endif
