// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SSHT_ADJOINT
#define SSHT_ADJOINT

#include <complex.h>

#ifdef __cplusplus
extern "C"{
#endif

void ssht_adjoint_mw_inverse_sov_sym(SSHT_COMPLEX(double) *flm, 
				     SSHT_COMPLEX(double) *f, 
				     int L, int spin, 
				     ssht_dl_method_t dl_method,
				     int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_real(SSHT_COMPLEX(double) *flm, 
					  double *f, 
					  int L,
					  ssht_dl_method_t dl_method, 
					  int verbosity);
void ssht_adjoint_mw_forward_sov_sym(SSHT_COMPLEX(double) *f, 
				     SSHT_COMPLEX(double) *flm,
				     int L, int spin,
				     ssht_dl_method_t dl_method,
				     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_real(double *f, 
					  SSHT_COMPLEX(double) *flm,
					  int L,
					  ssht_dl_method_t dl_method,
					  int verbosity);

void ssht_adjoint_mw_inverse_sov_sym_pole(SSHT_COMPLEX(double) *flm, SSHT_COMPLEX(double) *f,
					  SSHT_COMPLEX(double) f_sp, double phi_sp,
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_real_pole(SSHT_COMPLEX(double) *flm, 
					       double *f, 
					       double f_sp,
					       int L, 
					       ssht_dl_method_t dl_method,
					       int verbosity);
void ssht_adjoint_mw_forward_sov_sym_pole(SSHT_COMPLEX(double) *f, 
					  SSHT_COMPLEX(double) *f_sp, double *phi_sp,
					  SSHT_COMPLEX(double) *flm, 
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity);
void ssht_adjoint_mw_forward_sov_sym_real_pole(double *f, 
					       double *f_sp,
					       SSHT_COMPLEX(double) *flm, 
					       int L, 
					       ssht_dl_method_t dl_method, 
					       int verbosity);


void ssht_adjoint_mw_inverse_sov_sym_ss(SSHT_COMPLEX(double) *flm, SSHT_COMPLEX(double) *f, 
					int L, int spin, 
					ssht_dl_method_t dl_method,
					int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_ss_real(SSHT_COMPLEX(double) *flm, double *f, 
					     int L, 
					     ssht_dl_method_t dl_method, 
					     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss(SSHT_COMPLEX(double) *f, SSHT_COMPLEX(double) *flm,
					int L, int spin,
					ssht_dl_method_t dl_method,
					int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_real(double *f, 
					     SSHT_COMPLEX(double) *flm,
					     int L,
					     ssht_dl_method_t dl_method,
					     int verbosity);

void ssht_adjoint_mw_inverse_sov_sym_ss_pole(SSHT_COMPLEX(double) *flm, SSHT_COMPLEX(double) *f,
					     SSHT_COMPLEX(double) f_np, double phi_np,
					     SSHT_COMPLEX(double) f_sp, double phi_sp,
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity);
void ssht_adjoint_mw_inverse_sov_sym_ss_real_pole(SSHT_COMPLEX(double) *flm, 
						  double *f, 
						  double f_np,
						  double f_sp,
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_pole(SSHT_COMPLEX(double) *f, 
					     SSHT_COMPLEX(double) *f_np, double *phi_np,
					     SSHT_COMPLEX(double) *f_sp, double *phi_sp,
					     SSHT_COMPLEX(double) *flm, 
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity);
void ssht_adjoint_mw_forward_sov_sym_ss_real_pole(double *f, 
						  double *f_np,
						  double *f_sp,
						  SSHT_COMPLEX(double) *flm, 
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity);

#ifdef __cplusplus
}
#endif

#endif
