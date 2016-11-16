// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SSHT_CORE
#define SSHT_CORE

#include <complex.h>



void ssht_core_mw_inverse_sov_sym(SSHT_COMPLEX(double) *f, const SSHT_COMPLEX(double) *flm,
				  int L, int spin,
				  ssht_dl_method_t dl_method,
				  int verbosity);
void ssht_core_mw_lb_inverse_sov_sym(SSHT_COMPLEX(double) *f, const SSHT_COMPLEX(double) *flm,
				  int L0, int L, int spin,
				  ssht_dl_method_t dl_method,
				  int verbosity);
void ssht_core_mw_inverse_sov_sym_real(double *f, const SSHT_COMPLEX(double) *flm,
				       int L,
				       ssht_dl_method_t dl_method,
				       int verbosity);
void ssht_core_mw_lb_inverse_sov_sym_real(double *f, const SSHT_COMPLEX(double) *flm,
				       int L0, int L,
				       ssht_dl_method_t dl_method,
				       int verbosity);
void ssht_core_mw_forward_sov_conv_sym(SSHT_COMPLEX(double) *flm, const SSHT_COMPLEX(double) *f,
				       int L, int spin,
				       ssht_dl_method_t dl_method,
				       int verbosity);
void ssht_core_mw_lb_forward_sov_conv_sym(SSHT_COMPLEX(double) *flm, const SSHT_COMPLEX(double) *f,
				       int L0, int L, int spin,
				       ssht_dl_method_t dl_method,
				       int verbosity);
void ssht_core_mw_forward_sov_conv_sym_real(SSHT_COMPLEX(double) *flm, const double *f,
					    int L,
					    ssht_dl_method_t dl_method,
					    int verbosity);
void ssht_core_mw_lb_forward_sov_conv_sym_real(SSHT_COMPLEX(double) *flm, const double *f,
					    int L0, int L,
					    ssht_dl_method_t dl_method,
					    int verbosity);
void ssht_core_mw_inverse_sov_sym_pole(SSHT_COMPLEX(double) *f,
				       SSHT_COMPLEX(double) *f_sp, double *phi_sp,
				       const SSHT_COMPLEX(double) *flm,
				       int L, int spin,
				       ssht_dl_method_t dl_method,
				       int verbosity);
void ssht_core_mw_inverse_sov_sym_real_pole(double *f,
					    double *f_sp,
					    const SSHT_COMPLEX(double) *flm,
					    int L,
					    ssht_dl_method_t dl_method,
					    int verbosity);
void ssht_core_mw_forward_sov_conv_sym_pole(SSHT_COMPLEX(double) *flm, const SSHT_COMPLEX(double) *f,
					    SSHT_COMPLEX(double) f_sp, double phi_sp,
					    int L, int spin,
					    ssht_dl_method_t dl_method,
					    int verbosity);
void ssht_core_mw_forward_sov_conv_sym_real_pole(SSHT_COMPLEX(double) *flm,
						 const double *f,
						 double f_sp,
						 int L,
						 ssht_dl_method_t dl_method,
						 int verbosity);
// Note that mw direct algoritms are for testing purposes only.
void ssht_core_mwdirect_inverse(SSHT_COMPLEX(double) *f, const SSHT_COMPLEX(double) *flm,
				 int L, int spin, int verbosity);
void ssht_core_mwdirect_inverse_sov(SSHT_COMPLEX(double) *f, const SSHT_COMPLEX(double) *flm,
				     int L, int spin, int verbosity);


void ssht_core_mw_inverse_sov_sym_ss(SSHT_COMPLEX(double) *f, const SSHT_COMPLEX(double) *flm,
				     int L, int spin,
				     ssht_dl_method_t dl_method,
				     int verbosity);
void ssht_core_mw_lb_inverse_sov_sym_ss(SSHT_COMPLEX(double) *f, const SSHT_COMPLEX(double) *flm,
				     int L0, int L, int spin,
				     ssht_dl_method_t dl_method,
				     int verbosity);
void ssht_core_mw_inverse_sov_sym_ss_real(double *f, const SSHT_COMPLEX(double) *flm,
					  int L,
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_core_mw_lb_inverse_sov_sym_ss_real(double *f, const SSHT_COMPLEX(double) *flm,
					  int L0, int L,
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_core_mw_forward_sov_conv_sym_ss(SSHT_COMPLEX(double) *flm, const SSHT_COMPLEX(double) *f,
					  int L, int spin,
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_core_mw_lb_forward_sov_conv_sym_ss(SSHT_COMPLEX(double) *flm, const SSHT_COMPLEX(double) *f,
					  int L0, int L, int spin,
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_core_mw_forward_sov_conv_sym_ss_real(SSHT_COMPLEX(double) *flm, const double *f,
					       int L,
					       ssht_dl_method_t dl_method,
					       int verbosity);
void ssht_core_mw_lb_forward_sov_conv_sym_ss_real(SSHT_COMPLEX(double) *flm, const double *f,
					       int L0, int L,
					       ssht_dl_method_t dl_method,
					       int verbosity);
void ssht_core_mw_inverse_sov_sym_ss_pole(SSHT_COMPLEX(double) *f,
					  SSHT_COMPLEX(double) *f_np, double *phi_np,
					  SSHT_COMPLEX(double) *f_sp, double *phi_sp,
					  const SSHT_COMPLEX(double) *flm,
					  int L, int spin,
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_core_mw_inverse_sov_sym_ss_real_pole(double *f,
					       double *f_np,
					       double *f_sp,
					       const SSHT_COMPLEX(double) *flm,
					       int L,
					       ssht_dl_method_t dl_method,
					       int verbosity);
void ssht_core_mw_forward_sov_conv_sym_ss_pole(SSHT_COMPLEX(double) *flm, const SSHT_COMPLEX(double) *f,
					       SSHT_COMPLEX(double) f_np, double phi_np,
					       SSHT_COMPLEX(double) f_sp, double phi_sp,
					       int L, int spin,
					       ssht_dl_method_t dl_method,
					       int verbosity);
void ssht_core_mw_forward_sov_conv_sym_ss_real_pole(SSHT_COMPLEX(double) *flm,
						    const double *f,
						    double f_np,
						    double f_sp,
						    int L,
						    ssht_dl_method_t dl_method,
						    int verbosity);
// Note that mw direct algoritms are for testing purposes only.
void ssht_core_mwdirect_inverse_ss(SSHT_COMPLEX(double) *f, const SSHT_COMPLEX(double) *flm,
				   int L, int spin, int verbosity);


void ssht_core_gl_inverse_sov(SSHT_COMPLEX(double) *f, const SSHT_COMPLEX(double) *flm,
				     int L, int spin, int verbosity);
void ssht_core_gl_inverse_sov_real(double *f, const SSHT_COMPLEX(double) *flm,
				   int L, int verbosity);
void ssht_core_gl_forward_sov(SSHT_COMPLEX(double) *flm, const SSHT_COMPLEX(double) *f,
			      int L, int spin, int verbosity);
void ssht_core_gl_forward_sov_real(SSHT_COMPLEX(double) *flm, const double *f,
				   int L, int verbosity);


void ssht_core_dh_inverse_sov(SSHT_COMPLEX(double) *f, const SSHT_COMPLEX(double) *flm,
				     int L, int spin, int verbosity);
void ssht_core_dh_inverse_sov_real(double *f, const SSHT_COMPLEX(double) *flm,
				   int L, int verbosity);
void ssht_core_dh_forward_sov(SSHT_COMPLEX(double) *flm, const SSHT_COMPLEX(double) *f,
			      int L, int spin, int verbosity);
void ssht_core_dh_forward_sov_real(SSHT_COMPLEX(double) *flm, const double *f,
				   int L, int verbosity);


#endif
