// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SSHT_ADJOINT
#define SSHT_ADJOINT

#include "ssht_types.h"

#ifdef __cplusplus
extern "C"{
#endif

void ssht_adjoint_mw_inverse_sov_sym(ssht_complex_double *flm, 
				     const ssht_complex_double *f, 
				     int L, int spin, 
				     ssht_dl_method_t dl_method,
				     int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_real(ssht_complex_double *flm, 
					  const double *f, 
					  int L,
					  ssht_dl_method_t dl_method, 
					  int verbosity);
void ssht_adjoint_mw_forward_sov_sym(ssht_complex_double *f, 
				     const ssht_complex_double *flm,
				     int L, int spin,
				     ssht_dl_method_t dl_method,
				     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_real(double *f, 
					  const ssht_complex_double *flm,
					  int L,
					  ssht_dl_method_t dl_method,
					  int verbosity);

void ssht_adjoint_mw_inverse_sov_sym_pole(ssht_complex_double *flm, ssht_complex_double *f,
					  ssht_complex_double f_sp, double phi_sp,
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_real_pole(ssht_complex_double *flm, 
					       double *f, 
					       double f_sp,
					       int L, 
					       ssht_dl_method_t dl_method,
					       int verbosity);
void ssht_adjoint_mw_forward_sov_sym_pole(ssht_complex_double *f, 
					  ssht_complex_double *f_sp, double *phi_sp,
					  ssht_complex_double *flm, 
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_adjoint_mw_forward_sov_sym_real_pole(double *f, 
					       double *f_sp,
					       ssht_complex_double *flm, 
					       int L, 
					       ssht_dl_method_t dl_method, 
					       int verbosity);


void ssht_adjoint_mw_inverse_sov_sym_ss(ssht_complex_double *flm, ssht_complex_double *f, 
					int L, int spin, 
					ssht_dl_method_t dl_method,
					int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_ss_real(ssht_complex_double *flm, double *f, 
					     int L, 
					     ssht_dl_method_t dl_method, 
					     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss(ssht_complex_double *f, ssht_complex_double *flm,
					int L, int spin,
					ssht_dl_method_t dl_method,
					int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_real(double *f, 
					     ssht_complex_double *flm,
					     int L,
					     ssht_dl_method_t dl_method,
					     int verbosity);

void ssht_adjoint_mw_inverse_sov_sym_ss_pole(ssht_complex_double *flm, ssht_complex_double *f,
					     ssht_complex_double f_np, double phi_np,
					     ssht_complex_double f_sp, double phi_sp,
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_ss_real_pole(ssht_complex_double *flm, 
						  double *f, 
						  double f_np,
						  double f_sp,
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_pole(ssht_complex_double *f, 
					     ssht_complex_double *f_np, double *phi_np,
					     ssht_complex_double *f_sp, double *phi_sp,
					     ssht_complex_double *flm, 
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_real_pole(double *f, 
						  double *f_np,
						  double *f_sp,
						  ssht_complex_double *flm, 
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity);

#ifdef __cplusplus
}
#endif

#endif
