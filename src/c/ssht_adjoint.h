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
void ssht_adjoint_mw_inverse_sov_sym_real(complex double *flm, 
					  double *f, 
					  int L,
					  ssht_dl_method_t dl_method, 
					  int verbosity);
void ssht_adjoint_mw_forward_sov_sym(complex double *f, 
				     complex double *flm,
				     int L, int spin,
				     ssht_dl_method_t dl_method,
				     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_real(double *f, 
					  complex double *flm,
					  int L,
					  ssht_dl_method_t dl_method,
					  int verbosity);

void ssht_adjoint_mw_inverse_sov_sym_pole(complex double *flm, complex double *f,
					  complex double f_sp, double phi_sp,
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_real_pole(complex double *flm, 
					       double *f, 
					       double f_sp,
					       int L, 
					       ssht_dl_method_t dl_method,
					       int verbosity);
void ssht_adjoint_mw_forward_sov_sym_pole(complex double *f, 
					  complex double *f_sp, double *phi_sp,
					  complex double *flm, 
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_adjoint_mw_forward_sov_sym_real_pole(double *f, 
					       double *f_sp,
					       complex double *flm, 
					       int L, 
					       ssht_dl_method_t dl_method, 
					       int verbosity);


void ssht_adjoint_mw_inverse_sov_sym_ss(complex double *flm, complex double *f, 
					int L, int spin, 
					ssht_dl_method_t dl_method,
					int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_ss_real(complex double *flm, double *f, 
					     int L, 
					     ssht_dl_method_t dl_method, 
					     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss(complex double *f, complex double *flm,
					int L, int spin,
					ssht_dl_method_t dl_method,
					int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_real(double *f, 
					     complex double *flm,
					     int L,
					     ssht_dl_method_t dl_method,
					     int verbosity);

void ssht_adjoint_mw_inverse_sov_sym_ss_pole(complex double *flm, complex double *f,
					     complex double f_np, double phi_np,
					     complex double f_sp, double phi_sp,
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_ss_real_pole(complex double *flm, 
						  double *f, 
						  double f_np,
						  double f_sp,
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_pole(complex double *f, 
					     complex double *f_np, double *phi_np,
					     complex double *f_sp, double *phi_sp,
					     complex double *flm, 
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_real_pole(double *f, 
						  double *f_np,
						  double *f_sp,
						  complex double *flm, 
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity);


#endif
