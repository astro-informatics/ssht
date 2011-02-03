
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
// Note that mw direct algoritms are for testing purposes only.
void ssht_core_mwdirect_inverse(complex double *f, complex double *flm, 
				 int L, int spin, int verbosity);
void ssht_core_mwdirect_inverse_sov(complex double *f, complex double *flm, 
				     int L, int spin, int verbosity);


void ssht_core_mw_inverse_sov_sym_ss(complex double *f, complex double *flm, 
				     int L, int spin, int verbosity);
void ssht_core_mw_inverse_sov_sym_ss_real(double *f, complex double *flm, 
					  int L, int verbosity);
void ssht_core_mw_forward_sov_conv_sym_ss(complex double *flm, complex double *f, 
					  int L, int spin, int verbosity);







// Note that mw direct algoritms are for testing purposes only.
void ssht_core_mwdirect_inverse_ss(complex double *f, complex double *flm, 
				   int L, int spin, int verbosity);


void ssht_core_gl_inverse_sov(complex double *f, complex double *flm, 
				     int L, int spin, int verbosity);
void ssht_core_gl_inverse_sov_real(double *f, complex double *flm, 
				   int L, int verbosity);
void ssht_core_gl_forward_sov(complex double *flm, complex double *f, 
			      int L, int spin, int verbosity);
void ssht_core_gl_forward_sov_real(complex double *flm, double *f, 
				   int L, int verbosity);


void ssht_core_dh_inverse_sov(complex double *f, complex double *flm, 
				     int L, int spin, int verbosity);
void ssht_core_dh_inverse_sov_real(double *f, complex double *flm, 
				   int L, int verbosity);
void ssht_core_dh_forward_sov(complex double *flm, complex double *f, 
			      int L, int spin, int verbosity);
void ssht_core_dh_forward_sov_real(complex double *flm, double *f, 
				   int L, int verbosity);


#endif
