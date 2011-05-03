// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SSHT_ADJOINT
#define SSHT_ADJOINT

#include <complex.h>


void ssht_adjoint_mw_inverse_sov_sym(complex double *flm, 
				     complex double *f, 
				     int L, int spin, 
				     ssht_dl_method_t dl_method,
				     int verbosity);

void ssht_adjoint_mw_inverse_sov_sym_real(complex double *flm, double *f, 
					  int L,
					  ssht_dl_method_t dl_method, 
					  int verbosity);

#endif
