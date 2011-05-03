// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


/*! 
 * \file ssht_adjoint.c
 * Algorithms to perform fast adjoint of spin spherical harmonic
 * transforms on the sphere.
 *
 * \author Jason McEwen
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>  // Must be before fftw3.h
#include <fftw3.h>

#include "ssht_types.h"
#include "ssht_error.h"
#include "ssht_dl.h"
#include "ssht_sampling.h"
#include "ssht_adjoint.h"


//============================================================================
// MW algorithms
//============================================================================


/*!
 * Compute adjoint of inverse transform for MW method using separation
 * of variables, fast Fourier transforms and exploiting all symmetries
 * (for complex spin signal).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_inverse_sov_sym(complex double *flm, 
				     complex double *f, 
				     int L, int spin, 
				     ssht_dl_method_t dl_method,
				     int verbosity) {

  int el, m, mm, ind, t, r;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan, plan_bwd, plan_fwd;
  complex double *inout;
  complex double *Fmt, *Fmm, *Gmm;
  complex double *w, *wr;
  complex double *Fmm_pad, *tmp_pad;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  int w_offset;
  complex double *expsm, *expsmm;
  int exps_offset;
  int elmmsign, elssign;
  int spinneg;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  expsm = (complex double*)calloc(2*L-1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsm)
  expsmm = (complex double*)calloc(2*L-1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsmm)
  inds = (int*)calloc(2*L-1, sizeof(int));
  SSHT_ERROR_MEM_ALLOC_CHECK(inds)

  // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  for (m=0; m<=L-1; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }
  ssign = signs[abs(spin)];
  spinneg = spin <= 0 ? spin : -spin;
  exps_offset = L-1;
  for (m=-(L-1); m<=L-1; m++)
    expsm[m + exps_offset] = cexp(I*SSHT_PION2*(m+spin));
  for (mm=-(L-1); mm<=L-1; mm++)
    expsmm[mm + exps_offset] = cexp(-I*mm*SSHT_PI/(2.0*L-1.0));

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing forward transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_mw_forward_sov_conv_sym...");
  }

  // Compute Fourier transform over phi, i.e. compute Fmt.
  Fmt = (complex double*)calloc((2*L-1)*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmt)
  Fmt_stride = 2*L-1;
  Fmt_offset = L-1;
  f_stride = 2*L-1;
  inout = (complex double*)calloc(2*L-1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (t=0; t<=L-1; t++) {
    memcpy(inout, &f[t*f_stride], f_stride*sizeof(double complex));
    fftw_execute_dft(plan, inout, inout);
    for(m=0; m<=L-1; m++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m];
// / (2.0*L-1.0);
    for(m=-(L-1); m<=-1; m++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m+2*L-1]; 
// / (2.0*L-1.0);
  }

  // Extend Fmt periodically.
  for (m=-(L-1); m<=L-1; m++) 
    for (t=L; t<=2*L-2; t++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = 
	0.0;
//	signs[abs(m)] * ssign * Fmt[(m+Fmt_offset)*Fmt_stride + (2*L-2-t)];

  // Compute Fourier transform over theta, i.e. compute Fmm.
  Fmm = (complex double*)calloc((2*L-1)*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L-1;
  Fmm_offset = L-1;
  for (m=-(L-1); m<=L-1; m++) {
    memcpy(inout, &Fmt[(m+Fmt_offset)*Fmt_stride], Fmt_stride*sizeof(complex double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L-1; mm++) 
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] = 
	inout[mm]; 
// / (2.0*L-1.0);
    for(mm=-(L-1); mm<=-1; mm++) 
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] = 
	inout[mm+2*L-1]; 
// / (2.0*L-1.0);
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (m=-(L-1); m<=L-1; m++)
    for (mm=-(L-1); mm<=L-1; mm++)
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] *=
  	expsmm[mm + exps_offset];

  // Compute weights.
  w = (double complex*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(w)
  w_offset = 2*(L-1);
  for (mm=-2*(L-1); mm<=2*(L-1); mm++)
    w[mm+w_offset] = ssht_sampling_weight_mw(mm);

  // Compute IFFT of w to give wr.
  wr = (double complex*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(wr)
  inout = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan_bwd = fftw_plan_dft_1d(4*L-3, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  plan_fwd = fftw_plan_dft_1d(4*L-3, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (mm=1; mm<=2*L-2; mm++) 
    inout[mm + w_offset] = w[mm - 2*(L-1) - 1 + w_offset];
  for (mm=-2*(L-1); mm<=0; mm++) 
    inout[mm + w_offset] = w[mm + 2*(L-1) + w_offset];
  fftw_execute_dft(plan_bwd, inout, inout);
  for (mm=0; mm<=2*L-2; mm++) 
    wr[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
  for (mm=-2*(L-1); mm<=-1; mm++) 
    wr[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

  // Compute Gmm by convolution implemented as product in real space.
  Fmm_pad = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_pad)
  tmp_pad = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(tmp_pad)
  Gmm = (complex double*)calloc((2*L-1)*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Gmm)
  /* for (m=-(L-1); m<=L-1; m++) { */

  /*   // Zero-pad Fmm. */
  /*   for (mm=-2*(L-1); mm<=-L; mm++) */
  /*     Fmm_pad[mm+w_offset] = 0.0; */
  /*   for (mm=L; mm<=2*(L-1); mm++) */
  /*     Fmm_pad[mm+w_offset] = 0.0; */
  /*   for (mm=-(L-1); mm<=L-1; mm++) */
  /*     Fmm_pad[mm + w_offset] =  */
  /* 	Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset]; */

  /*   // Compute IFFT of Fmm. */
  /*   for (mm=1; mm<=2*L-2; mm++) */
  /*     inout[mm + w_offset] = Fmm_pad[mm - 2*(L-1) - 1 + w_offset]; */
  /*   for (mm=-2*(L-1); mm<=0; mm++) */
  /*     inout[mm + w_offset] = Fmm_pad[mm + 2*(L-1) + w_offset]; */
  /*   fftw_execute_dft(plan_bwd, inout, inout); */
  /*   for (mm=0; mm<=2*L-2; mm++) */
  /*     Fmm_pad[mm + w_offset] = inout[mm - 2*(L-1) + w_offset]; */
  /*   for (mm=-2*(L-1); mm<=-1; mm++) */
  /*     Fmm_pad[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset]; */

  /*   // Compute product of Fmm and weight in real space. */
  /*   for (r=-2*(L-1); r<=2*(L-1); r++)  */
  /*     Fmm_pad[r + w_offset] *= wr[-r + w_offset]; */

  /*   // Compute Gmm by FFT. */
  /*   for (mm=1; mm<=2*L-2; mm++) */
  /*     inout[mm + w_offset] = Fmm_pad[mm - 2*(L-1) - 1 + w_offset]; */
  /*   for (mm=-2*(L-1); mm<=0; mm++) */
  /*     inout[mm + w_offset] = Fmm_pad[mm + 2*(L-1) + w_offset]; */
  /*   fftw_execute_dft(plan_fwd, inout, inout); */
  /*   for (mm=0; mm<=2*L-2; mm++) */
  /*     Fmm_pad[mm + w_offset] = inout[mm - 2*(L-1) + w_offset]; */
  /*   for (mm=-2*(L-1); mm<=-1; mm++) */
  /*     Fmm_pad[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset]; */

  /*   // Extract section of Gmm of interest. */
  /*   for (mm=-(L-1); mm<=L-1; mm++) */
  /*     Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] =  */
  /* 	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0); */

  /* } */
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);

  // Compute flm.
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);   
  inds_offset = L-1;
  for (el=0; el<=L-1; el++) {
    for (m=-el; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (el=abs(spin); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==abs(spin)) {
	  for(eltmp=0; eltmp<=abs(spin); eltmp++)
	    ssht_dl_beta_risbo_eighth_table(dl8, SSHT_PION2, L, 
					    SSHT_DL_QUARTER_EXTENDED,
					    eltmp, sqrt_tbl, signs);
	  ssht_dl_beta_risbo_fill_eighth2quarter_table(dl, 
						       dl8, L,
						       SSHT_DL_QUARTER,
						       SSHT_DL_QUARTER_EXTENDED,
						       el, 
						       signs);
	}
	else {
	  ssht_dl_beta_risbo_eighth_table(dl8, SSHT_PION2, L, 
					  SSHT_DL_QUARTER_EXTENDED,
					  el, sqrt_tbl, signs);
	  ssht_dl_beta_risbo_fill_eighth2quarter_table(dl, 
						       dl8, L,
						       SSHT_DL_QUARTER,
						       SSHT_DL_QUARTER_EXTENDED,
						       el, 
						       signs);
	}
	break;
  
      case SSHT_DL_TRAPANI:
	if (el!=0 && el==abs(spin)) {
	  for(eltmp=0; eltmp<=abs(spin); eltmp++)
	    ssht_dl_halfpi_trapani_eighth_table(dl, L,
						SSHT_DL_QUARTER,
						eltmp, sqrt_tbl);
	  ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
							   SSHT_DL_QUARTER,
							   el, signs);
	}
	else {
	  ssht_dl_halfpi_trapani_eighth_table(dl, L,
					      SSHT_DL_QUARTER,
					      el, sqrt_tbl);
	  ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
							   SSHT_DL_QUARTER,
							   el, signs);	
	}
	break;

      default:
	SSHT_ERROR_GENERIC("Invalid dl method") 
    }

    // Compute flm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el *el + el;    
    for (m=-el; m<=el; m++)
      inds[m + inds_offset] = el2pel + m; 
    elssign = spin <= 0 ? 1.0 : signs[el];

    for (m=-el; m<=-1; m++) {
      // mm = 0
      ind = inds[m + inds_offset];
      flm[ind] +=
	ssign 
	* elfactor
	* expsm[m + exps_offset]
	* signs[el] * dl[0*dl_stride - m + dl_offset]
	* elssign * dl[0*dl_stride - spinneg + dl_offset]
	* Fmm[(m+Fmm_offset)*Fmm_stride + 0 + Fmm_offset];
//	* Gmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
    }
    for (m=0; m<=el; m++) {
      // mm = 0
      ind = inds[m + inds_offset];
      flm[ind] +=
	ssign 
	* elfactor
	* expsm[m + exps_offset]
	* dl[0*dl_stride + m + dl_offset]
	* elssign * dl[0*dl_stride - spinneg + dl_offset]
	* Fmm[(m+Fmm_offset)*Fmm_stride + 0 + Fmm_offset];
//      * Gmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
    }

    for (mm=1; mm<=el; mm++) {
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      for (m=-el; m<=-1; m++) {
      	ind = inds[m + inds_offset];
      	flm[ind] +=
      	  ssign
      	  * elfactor
      	  * expsm[m + exps_offset]
      	  * elmmsign * dl[mm*dl_stride - m + dl_offset]
      	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
      	  * ( Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset]
      	      + signs[-m] * ssign
      	      * Fmm[(m+Fmm_offset)*Fmm_stride - mm + Fmm_offset]);
	  /* * ( Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] */
      	  /*     + signs[-m] * ssign */
      	  /*     * Gmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]); */
      }
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign 
	  * elfactor
	  * expsm[m + exps_offset]
	  * dl[mm*dl_stride + m + dl_offset]
	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
	  * ( Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset]
	      + signs[m] * ssign
	      * Fmm[(m+Fmm_offset)*Fmm_stride - mm + Fmm_offset]);
	/* * ( Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] */
	/*       + signs[m] * ssign */
	/*       * Gmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]); */
      }

    }  

  }

  // Free memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);
  free(Fmt);
  free(Fmm);
  free(inout);
  free(w);
  free(wr);
  free(Fmm_pad);
  free(tmp_pad);
  free(Gmm);
  free(sqrt_tbl);
  free(signs); 
  free(expsm);
  free(expsmm);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s %s", SSHT_PROMPT, "Forward transform computed!");  

}


/*!  
 * Compute adjoint of inverse transform for MW method of real scalar
 * signal using separation of variables, fast Fourier transforms and
 * exploiting all symmetries (including additional symmetries for real
 * signals).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_inverse_sov_sym_real(complex double *flm, double *f, 
					  int L,
					  ssht_dl_method_t dl_method, 
					  int verbosity) {

  int el, m, mm, ind, ind_nm, t, r;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan, plan_bwd, plan_fwd;
  double *in_real;
  complex double *inout, *out;
  complex double *Fmt, *Fmm, *Gmm;
  complex double *w, *wr;
  complex double *Fmm_pad, *tmp_pad;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset, Gmm_stride;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  int w_offset;
  complex double *expsm, *expsmm;
  int exps_offset;
  int elmmsign, elssign;
  int spinneg;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  expsm = (complex double*)calloc(L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsm)
  expsmm = (complex double*)calloc(2*L-1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsmm)
  inds = (int*)calloc(L, sizeof(int));
  SSHT_ERROR_MEM_ALLOC_CHECK(inds)

  // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  for (m=0; m<=L-1; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }
  ssign = signs[abs(spin)];
  spinneg = spin <= 0 ? spin : -spin;
  exps_offset = L-1;
  for (m=0; m<=L-1; m++)
    expsm[m] = cexp(I*SSHT_PION2*(m+spin));
  for (mm=-(L-1); mm<=L-1; mm++)
    expsmm[mm + exps_offset] = cexp(-I*mm*SSHT_PI/(2.0*L-1.0));

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing forward transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_mw_forward_sov_conv_sym_real...");
  }

  // Compute Fourier transform over phi, i.e. compute Fmt.
  Fmt = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmt)
  Fmt_stride = 2*L-1;
  Fmt_offset = L-1;
  f_stride = 2*L-1;
  in_real = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in_real)
  out = (complex double*)calloc(L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out)
  plan = fftw_plan_dft_r2c_1d(2*L-1, in_real, out, FFTW_MEASURE);
  for (t=0; t<=L-1; t++) {
    memcpy(in_real, &f[t*f_stride], f_stride*sizeof(double));
        fftw_execute_dft_r2c(plan, in_real, out);
    for(m=0; m<=L-1; m++) 
      Fmt[m*Fmt_stride + t] = out[m];
// / (2.0*L-1.0);
  }
  free(in_real);
  free(out);
  fftw_destroy_plan(plan);

  // Extend Fmt periodically.
  for (m=0; m<=L-1; m++) 
    for (t=L; t<=2*L-2; t++) 
      Fmt[m*Fmt_stride + t] = 
	0.0;
//	signs[abs(m)] * ssign * Fmt[m*Fmt_stride + (2*L-2-t)];

  // Compute Fourier transform over theta, i.e. compute Fmm.
  Fmm = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L-1;
  Fmm_offset = L-1;
  inout = (complex double*)calloc(2*L-1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (m=0; m<=L-1; m++) {
    memcpy(inout, &Fmt[m*Fmt_stride], Fmt_stride*sizeof(complex double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L-1; mm++) 
      Fmm[m*Fmm_stride + mm + Fmm_offset] = 
	inout[mm]; 
// / (2.0*L-1.0);
    for(mm=-(L-1); mm<=-1; mm++) 
      Fmm[m*Fmm_stride + mm + Fmm_offset] = 
	inout[mm+2*L-1]; 
// / (2.0*L-1.0);
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (m=0; m<=L-1; m++)
    for (mm=-(L-1); mm<=L-1; mm++)
      Fmm[m*Fmm_stride + mm + Fmm_offset] *= 
	expsmm[mm + exps_offset];

  // Compute weights.
  w = (double complex*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(w)
  w_offset = 2*(L-1);
  for (mm=-2*(L-1); mm<=2*(L-1); mm++)
    w[mm+w_offset] = ssht_sampling_weight_mw(mm);

  // Compute IFFT of w to give wr.
  wr = (double complex*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(wr)
  inout = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan_bwd = fftw_plan_dft_1d(4*L-3, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  plan_fwd = fftw_plan_dft_1d(4*L-3, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (mm=1; mm<=2*L-2; mm++) 
    inout[mm + w_offset] = w[mm - 2*(L-1) - 1 + w_offset];
  for (mm=-2*(L-1); mm<=0; mm++) 
    inout[mm + w_offset] = w[mm + 2*(L-1) + w_offset];
  fftw_execute_dft(plan_bwd, inout, inout);
  for (mm=0; mm<=2*L-2; mm++) 
    wr[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
  for (mm=-2*(L-1); mm<=-1; mm++) 
    wr[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

  // Compute Gmm by convolution implemented as product in real space.
  Fmm_pad = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_pad)
  tmp_pad = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(tmp_pad)
  Gmm = (complex double*)calloc((2*L-1)*L, sizeof(complex double));
  Gmm_stride = L;
  SSHT_ERROR_MEM_ALLOC_CHECK(Gmm)
  for (m=0; m<=L-1; m++) {

    // Zero-pad Fmm.
    for (mm=-2*(L-1); mm<=-L; mm++)
      Fmm_pad[mm+w_offset] = 0.0;
    for (mm=L; mm<=2*(L-1); mm++)
      Fmm_pad[mm+w_offset] = 0.0;
    for (mm=-(L-1); mm<=L-1; mm++)
      Fmm_pad[mm + w_offset] = 
	Fmm[m*Fmm_stride + mm + Fmm_offset];

    // Compute IFFT of Fmm.
    for (mm=1; mm<=2*L-2; mm++)
      inout[mm + w_offset] = Fmm_pad[mm - 2*(L-1) - 1 + w_offset];
    for (mm=-2*(L-1); mm<=0; mm++)
      inout[mm + w_offset] = Fmm_pad[mm + 2*(L-1) + w_offset];
    fftw_execute_dft(plan_bwd, inout, inout);
    for (mm=0; mm<=2*L-2; mm++)
      Fmm_pad[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
    for (mm=-2*(L-1); mm<=-1; mm++)
      Fmm_pad[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

    // Compute product of Fmm and weight in real space.
    for (r=-2*(L-1); r<=2*(L-1); r++) 
      Fmm_pad[r + w_offset] *= wr[-r + w_offset];

    // Compute Gmm by FFT.
    for (mm=1; mm<=2*L-2; mm++)
      inout[mm + w_offset] = Fmm_pad[mm - 2*(L-1) - 1 + w_offset];
    for (mm=-2*(L-1); mm<=0; mm++)
      inout[mm + w_offset] = Fmm_pad[mm + 2*(L-1) + w_offset];
    fftw_execute_dft(plan_fwd, inout, inout);
    for (mm=0; mm<=2*L-2; mm++)
      Fmm_pad[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
    for (mm=-2*(L-1); mm<=-1; mm++)
      Fmm_pad[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

    // Extract section of Gmm of interest.
    for (mm=-(L-1); mm<=L-1; mm++)
      Gmm[(mm+Fmm_offset)*Gmm_stride + m] = 
	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);

  }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);

  // Compute flm.
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER); 
  inds_offset = 0;
  for (el=0; el<=L-1; el++) {
    for (m=0; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (el=abs(spin); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==abs(spin)) {
	  for(eltmp=0; eltmp<=abs(spin); eltmp++)
	    ssht_dl_beta_risbo_eighth_table(dl8, SSHT_PION2, L, 
					    SSHT_DL_QUARTER_EXTENDED,
					    eltmp, sqrt_tbl, signs);
	  ssht_dl_beta_risbo_fill_eighth2quarter_table(dl, 
						       dl8, L,
						       SSHT_DL_QUARTER,
						       SSHT_DL_QUARTER_EXTENDED,
						       el, 
						       signs);
	}
	else {
	  ssht_dl_beta_risbo_eighth_table(dl8, SSHT_PION2, L, 
					  SSHT_DL_QUARTER_EXTENDED,
					  el, sqrt_tbl, signs);
	  ssht_dl_beta_risbo_fill_eighth2quarter_table(dl, 
						       dl8, L,
						       SSHT_DL_QUARTER,
						       SSHT_DL_QUARTER_EXTENDED,
						       el, 
						       signs);
	}
	break;
  
      case SSHT_DL_TRAPANI:
	if (el!=0 && el==abs(spin)) {
	  for(eltmp=0; eltmp<=abs(spin); eltmp++)
	    ssht_dl_halfpi_trapani_eighth_table(dl, L,
						SSHT_DL_QUARTER,
						eltmp, sqrt_tbl);
	  ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
							   SSHT_DL_QUARTER,
							   el, signs);
	}
	else {
	  ssht_dl_halfpi_trapani_eighth_table(dl, L,
					      SSHT_DL_QUARTER,
					      el, sqrt_tbl);
	  ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
							   SSHT_DL_QUARTER,
							   el, signs);	
	}
	break;

      default:
	SSHT_ERROR_GENERIC("Invalid dl method") 
    }


    // Compute flm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el *el + el;    
    for (m=0; m<=el; m++)
      inds[m + inds_offset] = el2pel + m; 
    elssign = spin <= 0 ? 1.0 : signs[el];

    for (m=0; m<=el; m++) {
      // mm = 0
      ind = inds[m + inds_offset];
      flm[ind] +=
	ssign 
	* elfactor
	* expsm[m]
	* dl[0*dl_stride + m + dl_offset]
	* elssign * dl[0*dl_stride - spinneg + dl_offset]
	* Fmm[m*Fmm_stride + 0 + Fmm_offset];
//	* Gmm[(0+Fmm_offset)*Gmm_stride + m];
    }

    for (mm=1; mm<=el; mm++) {
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign 
	  * elfactor
	  * expsm[m]
	  * dl[mm*dl_stride + m + dl_offset]
	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
	  * (  Fmm[m*Fmm_stride + mm + Fmm_offset] 
	      + signs[m] * ssign
	       * Fmm[m*Fmm_stride - mm + Fmm_offset]);
  	  /* * ( Gmm[(mm+Fmm_offset)*Gmm_stride + m] */
	  /*     + signs[m] * ssign */
	  /*     * Gmm[(-mm+Fmm_offset)*Gmm_stride + m]); */
      }

    }  

  }

  // Set flm values for negative m using conjugate symmetry.
  for (el=abs(spin); el<=L-1; el++) {
    for (m=1; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      ssht_sampling_elm2ind(&ind_nm, el, -m);
      flm[ind_nm] = signs[m] * conj(flm[ind]);
    }
  }

  // Free memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);
  free(Fmt);
  free(Fmm);
  free(inout);
  free(w);
  free(wr);
  free(Fmm_pad);
  free(tmp_pad);
  free(Gmm);
  free(sqrt_tbl);
  free(signs); 
  free(expsm);
  free(expsmm);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s %s", SSHT_PROMPT, "Forward transform computed!");  

}


//============================================================================
// MW SS algorithms
//============================================================================


/*!  
 * Compute adjoint of inverse transform for MW method with symmetric
 * sampling using separation of variables, fast Fourier transforms and
 * exploiting all symmetries (for complex spin signal).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_inverse_sov_sym_ss(complex double *flm, complex double *f, 
					int L, int spin, 
					ssht_dl_method_t dl_method,
					int verbosity) {

  int el, m, mm, ind, t, r;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan, plan_bwd, plan_fwd;
  complex double *inout;
  complex double *Fmt, *Fmm, *Gmm;
  complex double *w, *wr;
  complex double *Fmm_pad, *tmp_pad;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  int w_offset;
  complex double *expsm, *expsmm;
  int exps_offset;
  int elmmsign, elssign;
  int spinneg;
  int Gmm_stride, Gmm_offset;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  expsm = (complex double*)calloc(2*L-1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsm)
  expsmm = (complex double*)calloc(2*L-1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsmm)
  inds = (int*)calloc(2*L-1, sizeof(int));
  SSHT_ERROR_MEM_ALLOC_CHECK(inds)

  // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  for (m=0; m<=L-1; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }
  ssign = signs[abs(spin)];
  spinneg = spin <= 0 ? spin : -spin;
  exps_offset = L-1;
  for (m=-(L-1); m<=L-1; m++)
    expsm[m + exps_offset] = cexp(I*SSHT_PION2*(m+spin));
  for (mm=-(L-1); mm<=L-1; mm++)
    expsmm[mm + exps_offset] = cexp(-I*mm*SSHT_PI/(2.0*L-1.0));

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing forward transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_mw_forward_sov_conv_sym_ss...");
  }

  // Compute Fourier transform over phi, i.e. compute Fmt.
  // Note that t and p indices of fext are increased in size by
  // one compared to usual sampling.
  Fmt = (complex double*)calloc((2*L)*(2*L), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmt)
  Fmt_stride = 2*L;
  Fmt_offset = L-1;
  f_stride = 2*L;
  inout = (complex double*)calloc(2*L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (t=0; t<=L; t++) {
    memcpy(inout, &f[t*f_stride], f_stride*sizeof(double complex));
    fftw_execute_dft(plan, inout, inout);
    for(m=0; m<=L; m++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m]; 
// / (2.0*L);
    for(m=-(L-1); m<=-1; m++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m+2*L-1+1]; 
// / (2.0*L);
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Extend Fmt periodically.
  for (m=-(L-1); m<=L; m++) 
    for (t=L+1; t<=2*L-1; t++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = 
	0.0;
//	signs[abs(m)] * ssign * Fmt[(m+Fmt_offset)*Fmt_stride + (2*L-t)];

  // Compute Fourier transform over theta, i.e. compute Fmm.
  // Note that m and mm indices are increased in size by one.
  Fmm = (complex double*)calloc((2*L)*(2*L), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L;
  Fmm_offset = L-1;
  inout = (complex double*)calloc(2*L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (m=-(L-1); m<=L; m++) {
    memcpy(inout, &Fmt[(m+Fmt_offset)*Fmt_stride], Fmt_stride*sizeof(complex double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L; mm++) 
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] = 
	inout[mm]; 
// / (2.0*L);
    for(mm=-(L-1); mm<=-1; mm++) 
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] = 
	inout[mm+2*L-1+1]; 
// / (2.0*L);
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Compute weights.
  w = (double complex*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(w)
  w_offset = 2*(L-1);
  for (mm=-2*(L-1); mm<=2*(L-1); mm++)
    w[mm+w_offset] = ssht_sampling_weight_mw(mm);

  // Compute IFFT of w to give wr.
  wr = (double complex*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(wr)
  inout = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan_bwd = fftw_plan_dft_1d(4*L-3, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  plan_fwd = fftw_plan_dft_1d(4*L-3, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (mm=1; mm<=2*L-2; mm++) 
    inout[mm + w_offset] = w[mm - 2*(L-1) - 1 + w_offset];
  for (mm=-2*(L-1); mm<=0; mm++) 
    inout[mm + w_offset] = w[mm + 2*(L-1) + w_offset];
  fftw_execute_dft(plan_bwd, inout, inout);
  for (mm=0; mm<=2*L-2; mm++) 
    wr[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
  for (mm=-2*(L-1); mm<=-1; mm++) 
    wr[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

  // Compute Gmm by convolution implemented as product in real space.
  Fmm_pad = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_pad)
  tmp_pad = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(tmp_pad)
  Gmm = (complex double*)calloc((2*L-1)*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Gmm)
  Gmm_stride = 2*L-1;
  Gmm_offset = L-1;
  for (m=-(L-1); m<=L-1; m++) {

    // Zero-pad Fmm.
    for (mm=-2*(L-1); mm<=-L; mm++)
      Fmm_pad[mm+w_offset] = 0.0;
    for (mm=L; mm<=2*(L-1); mm++)
      Fmm_pad[mm+w_offset] = 0.0;
    for (mm=-(L-1); mm<=L-1; mm++)
      Fmm_pad[mm + w_offset] = 
	Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset];

    // Compute IFFT of Fmm.
    for (mm=1; mm<=2*L-2; mm++)
      inout[mm + w_offset] = Fmm_pad[mm - 2*(L-1) - 1 + w_offset];
    for (mm=-2*(L-1); mm<=0; mm++)
      inout[mm + w_offset] = Fmm_pad[mm + 2*(L-1) + w_offset];
    fftw_execute_dft(plan_bwd, inout, inout);
    for (mm=0; mm<=2*L-2; mm++)
      Fmm_pad[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
    for (mm=-2*(L-1); mm<=-1; mm++)
      Fmm_pad[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

    // Compute product of Fmm and weight in real space.
    for (r=-2*(L-1); r<=2*(L-1); r++) 
      Fmm_pad[r + w_offset] *= wr[-r + w_offset];

    // Compute Gmm by FFT.
    for (mm=1; mm<=2*L-2; mm++)
      inout[mm + w_offset] = Fmm_pad[mm - 2*(L-1) - 1 + w_offset];
    for (mm=-2*(L-1); mm<=0; mm++)
      inout[mm + w_offset] = Fmm_pad[mm + 2*(L-1) + w_offset];
    fftw_execute_dft(plan_fwd, inout, inout);
    for (mm=0; mm<=2*L-2; mm++)
      Fmm_pad[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
    for (mm=-2*(L-1); mm<=-1; mm++)
      Fmm_pad[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

    // Extract section of Gmm of interest.
    for (mm=-(L-1); mm<=L-1; mm++)
      Gmm[(mm+Gmm_offset)*Gmm_stride + m + Gmm_offset] = 
	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);

  }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);

  // Compute flm.
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER); 
  inds_offset = L-1;
  for (el=0; el<=L-1; el++) {
    for (m=-el; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (el=abs(spin); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==abs(spin)) {
	  for(eltmp=0; eltmp<=abs(spin); eltmp++)
	    ssht_dl_beta_risbo_eighth_table(dl8, SSHT_PION2, L, 
					    SSHT_DL_QUARTER_EXTENDED,
					    eltmp, sqrt_tbl, signs);
	  ssht_dl_beta_risbo_fill_eighth2quarter_table(dl, 
						       dl8, L,
						       SSHT_DL_QUARTER,
						       SSHT_DL_QUARTER_EXTENDED,
						       el, 
						       signs);
	}
	else {
	  ssht_dl_beta_risbo_eighth_table(dl8, SSHT_PION2, L, 
					  SSHT_DL_QUARTER_EXTENDED,
					  el, sqrt_tbl, signs);
	  ssht_dl_beta_risbo_fill_eighth2quarter_table(dl, 
						       dl8, L,
						       SSHT_DL_QUARTER,
						       SSHT_DL_QUARTER_EXTENDED,
						       el, 
						       signs);
	}
	break;
  
      case SSHT_DL_TRAPANI:
	if (el!=0 && el==abs(spin)) {
	  for(eltmp=0; eltmp<=abs(spin); eltmp++)
	    ssht_dl_halfpi_trapani_eighth_table(dl, L,
						SSHT_DL_QUARTER,
						eltmp, sqrt_tbl);
	  ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
							   SSHT_DL_QUARTER,
							   el, signs);
	}
	else {
	  ssht_dl_halfpi_trapani_eighth_table(dl, L,
					      SSHT_DL_QUARTER,
					      el, sqrt_tbl);
	  ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
							   SSHT_DL_QUARTER,
							   el, signs);	
	}
	break;

      default:
	SSHT_ERROR_GENERIC("Invalid dl method") 
    }

    // Compute flm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el *el + el;    
    for (m=-el; m<=el; m++)
      inds[m + inds_offset] = el2pel + m; 
    elssign = spin <= 0 ? 1.0 : signs[el];

    for (m=-el; m<=-1; m++) {
      // mm = 0
      ind = inds[m + inds_offset];
      flm[ind] +=
	ssign 
	* elfactor
	* expsm[m + exps_offset]
	* signs[el] * dl[0*dl_stride - m + dl_offset]
	* elssign * dl[0*dl_stride - spinneg + dl_offset]
	* Fmm[(m+Fmm_offset)*Fmm_stride + 0 + Fmm_offset];
//	* Gmm[(0+Gmm_offset)*Gmm_stride + m + Gmm_offset];
    }
    for (m=0; m<=el; m++) {
      // mm = 0
      ind = inds[m + inds_offset];
      flm[ind] +=
	ssign 
	* elfactor
	* expsm[m + exps_offset]
	* dl[0*dl_stride + m + dl_offset]
	* elssign * dl[0*dl_stride - spinneg + dl_offset]
	* Fmm[(m+Fmm_offset)*Fmm_stride + 0 + Fmm_offset];
//	* Gmm[(0+Gmm_offset)*Gmm_stride + m + Gmm_offset];
    }

    for (mm=1; mm<=el; mm++) {
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      for (m=-el; m<=-1; m++) {
      	ind = inds[m + inds_offset];
      	flm[ind] +=
      	  ssign
      	  * elfactor
      	  * expsm[m + exps_offset]
      	  * elmmsign * dl[mm*dl_stride - m + dl_offset]
      	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
      	  * ( Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset]
      	      + signs[-m] * ssign
      	      * Fmm[(m+Fmm_offset)*Fmm_stride - mm + Fmm_offset]);
      	  /* * ( Gmm[(mm+Gmm_offset)*Gmm_stride + m + Gmm_offset] */
      	  /*     + signs[-m] * ssign */
      	  /*     * Gmm[(-mm+Gmm_offset)*Gmm_stride + m + Gmm_offset]); */
      }
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign 
	  * elfactor
	  * expsm[m + exps_offset]
	  * dl[mm*dl_stride + m + dl_offset]
	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
	  * ( Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset]
	      + signs[m] * ssign
	      * Fmm[(m+Fmm_offset)*Fmm_stride - mm + Fmm_offset]);
	  /* * ( Gmm[(mm+Gmm_offset)*Gmm_stride + m + Gmm_offset] */
	  /*     + signs[m] * ssign */
	  /*     * Gmm[(-mm+Gmm_offset)*Gmm_stride + m + Gmm_offset]); */
      }

    }  

  }

  // Free memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);
  free(Fmt);
  free(Fmm);
  free(inout);
  free(w);
  free(wr);
  free(Fmm_pad);
  free(tmp_pad);
  free(Gmm);
  free(sqrt_tbl);
  free(signs); 
  free(expsm);
  free(expsmm);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s %s", SSHT_PROMPT, "Forward transform computed!");  

}


/*!  
 * Compute adjoint of inverse transform for MW method with symmetric
 * sampling of real scalar signal using separation of variables, fast
 * Fourier transforms and exploiting all symmetries (including
 * additional symmetries for real signals).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_inverse_sov_sym_ss_real(complex double *flm, double *f, 
					     int L, 
					     ssht_dl_method_t dl_method, 
					     int verbosity) {

  int el, m, mm, ind, ind_nm, t, r;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan, plan_bwd, plan_fwd;
  double *in_real;
  complex double *inout, *out;
  complex double *Fmt, *Fmm, *Gmm;
  complex double *w, *wr;
  complex double *Fmm_pad, *tmp_pad;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset, Gmm_stride;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  int w_offset;
  complex double *expsm;
  int exps_offset;
  int elmmsign, elssign;
  int spinneg;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  expsm = (complex double*)calloc(L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsm)
  inds = (int*)calloc(L, sizeof(int));
  SSHT_ERROR_MEM_ALLOC_CHECK(inds)

  // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  for (m=0; m<=L-1; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }
  ssign = signs[abs(spin)];
  spinneg = spin <= 0 ? spin : -spin;
  exps_offset = L-1;
  for (m=0; m<=L-1; m++)
    expsm[m] = cexp(I*SSHT_PION2*(m+spin));

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing forward transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_mw_forward_sov_conv_sym_ss_real...");
  }

  // Compute Fourier transform over phi, i.e. compute Fmt.
  // Note that t and p indices of fext are increased in size by
  // one compared to usual sampling.
  Fmt = (complex double*)calloc((L+1)*(2*L), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmt)
  Fmt_stride = 2*L;
  Fmt_offset = L-1;
  f_stride = 2*L;
  in_real = (double*)calloc(2*L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in_real)
  out = (complex double*)calloc(L+1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out)
  plan = fftw_plan_dft_r2c_1d(2*L, in_real, out, FFTW_MEASURE);
  for (t=0; t<=L; t++) {
    memcpy(in_real, &f[t*f_stride], f_stride*sizeof(double));
        fftw_execute_dft_r2c(plan, in_real, out);
    for(m=0; m<=L; m++) 
      Fmt[m*Fmt_stride + t] = out[m]; 
// / (2.0*L);
  }
  free(in_real);
  free(out);
  fftw_destroy_plan(plan);

  // Extend Fmt periodically.
  for (m=0; m<=L; m++) 
    for (t=L+1; t<=2*L-1; t++) 
      Fmt[m*Fmt_stride + t] = 
	0.0;
//	signs[abs(m)] * ssign * Fmt[m*Fmt_stride + (2*L-t)];

  // Compute Fourier transform over theta, i.e. compute Fmm.
  // Note that m and mm indices are increased in size by one.
  Fmm = (complex double*)calloc((L+1)*(2*L), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L;
  Fmm_offset = L-1;
  inout = (complex double*)calloc(2*L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (m=0; m<=L; m++) {
    memcpy(inout, &Fmt[m*Fmt_stride], Fmt_stride*sizeof(complex double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L; mm++) 
      Fmm[m*Fmm_stride + mm + Fmm_offset] = 
	inout[mm]; 
// / (2.0*L);
    for(mm=-(L-1); mm<=-1; mm++) 
      Fmm[m*Fmm_stride + mm + Fmm_offset] = 
	inout[mm+2*L-1+1]; 
// / (2.0*L);
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Compute weights.
  w = (double complex*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(w)
  w_offset = 2*(L-1);
  for (mm=-2*(L-1); mm<=2*(L-1); mm++)
    w[mm+w_offset] = ssht_sampling_weight_mw(mm);

  // Compute IFFT of w to give wr.
  wr = (double complex*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(wr)
  inout = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan_bwd = fftw_plan_dft_1d(4*L-3, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  plan_fwd = fftw_plan_dft_1d(4*L-3, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (mm=1; mm<=2*L-2; mm++) 
    inout[mm + w_offset] = w[mm - 2*(L-1) - 1 + w_offset];
  for (mm=-2*(L-1); mm<=0; mm++) 
    inout[mm + w_offset] = w[mm + 2*(L-1) + w_offset];
  fftw_execute_dft(plan_bwd, inout, inout);
  for (mm=0; mm<=2*L-2; mm++) 
    wr[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
  for (mm=-2*(L-1); mm<=-1; mm++) 
    wr[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

  // Compute Gmm by convolution implemented as product in real space.
  Fmm_pad = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_pad)
  tmp_pad = (complex double*)calloc(4*L-3, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(tmp_pad)
  Gmm = (complex double*)calloc((2*L-1)*L, sizeof(complex double));
  Gmm_stride = L;
  SSHT_ERROR_MEM_ALLOC_CHECK(Gmm)
  for (m=0; m<=L-1; m++) {

    // Zero-pad Fmm.
    for (mm=-2*(L-1); mm<=-L; mm++)
      Fmm_pad[mm+w_offset] = 0.0;
    for (mm=L; mm<=2*(L-1); mm++)
      Fmm_pad[mm+w_offset] = 0.0;
    for (mm=-(L-1); mm<=L-1; mm++)
      Fmm_pad[mm + w_offset] = 
	Fmm[m*Fmm_stride + mm + Fmm_offset];

    // Compute IFFT of Fmm.
    for (mm=1; mm<=2*L-2; mm++)
      inout[mm + w_offset] = Fmm_pad[mm - 2*(L-1) - 1 + w_offset];
    for (mm=-2*(L-1); mm<=0; mm++)
      inout[mm + w_offset] = Fmm_pad[mm + 2*(L-1) + w_offset];
    fftw_execute_dft(plan_bwd, inout, inout);
    for (mm=0; mm<=2*L-2; mm++)
      Fmm_pad[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
    for (mm=-2*(L-1); mm<=-1; mm++)
      Fmm_pad[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

    // Compute product of Fmm and weight in real space.
    for (r=-2*(L-1); r<=2*(L-1); r++) 
      Fmm_pad[r + w_offset] *= wr[-r + w_offset];

    // Compute Gmm by FFT.
    for (mm=1; mm<=2*L-2; mm++)
      inout[mm + w_offset] = Fmm_pad[mm - 2*(L-1) - 1 + w_offset];
    for (mm=-2*(L-1); mm<=0; mm++)
      inout[mm + w_offset] = Fmm_pad[mm + 2*(L-1) + w_offset];
    fftw_execute_dft(plan_fwd, inout, inout);
    for (mm=0; mm<=2*L-2; mm++)
      Fmm_pad[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
    for (mm=-2*(L-1); mm<=-1; mm++)
      Fmm_pad[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];

    // Extract section of Gmm of interest.
    for (mm=-(L-1); mm<=L-1; mm++)
      Gmm[(mm+Fmm_offset)*Gmm_stride + m] = 
	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);

  }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);

  // Compute flm.
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER); 
  inds_offset = 0;
  for (el=0; el<=L-1; el++) {
    for (m=0; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (el=abs(spin); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==abs(spin)) {
	  for(eltmp=0; eltmp<=abs(spin); eltmp++)
	    ssht_dl_beta_risbo_eighth_table(dl8, SSHT_PION2, L, 
					    SSHT_DL_QUARTER_EXTENDED,
					    eltmp, sqrt_tbl, signs);
	  ssht_dl_beta_risbo_fill_eighth2quarter_table(dl, 
						       dl8, L,
						       SSHT_DL_QUARTER,
						       SSHT_DL_QUARTER_EXTENDED,
						       el, 
						       signs);
	}
	else {
	  ssht_dl_beta_risbo_eighth_table(dl8, SSHT_PION2, L, 
					  SSHT_DL_QUARTER_EXTENDED,
					  el, sqrt_tbl, signs);
	  ssht_dl_beta_risbo_fill_eighth2quarter_table(dl, 
						       dl8, L,
						       SSHT_DL_QUARTER,
						       SSHT_DL_QUARTER_EXTENDED,
						       el, 
						       signs);
	}
	break;
  
      case SSHT_DL_TRAPANI:
	if (el!=0 && el==abs(spin)) {
	  for(eltmp=0; eltmp<=abs(spin); eltmp++)
	    ssht_dl_halfpi_trapani_eighth_table(dl, L,
						SSHT_DL_QUARTER,
						eltmp, sqrt_tbl);
	  ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
							   SSHT_DL_QUARTER,
							   el, signs);
	}
	else {
	  ssht_dl_halfpi_trapani_eighth_table(dl, L,
					      SSHT_DL_QUARTER,
					      el, sqrt_tbl);
	  ssht_dl_halfpi_trapani_fill_eighth2quarter_table(dl, L,
							   SSHT_DL_QUARTER,
							   el, signs);	
	}
	break;

      default:
	SSHT_ERROR_GENERIC("Invalid dl method") 
    }

    // Compute flm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el *el + el;    
    for (m=0; m<=el; m++)
      inds[m + inds_offset] = el2pel + m; 
    elssign = spin <= 0 ? 1.0 : signs[el];

    for (m=0; m<=el; m++) {
      // mm = 0
      ind = inds[m + inds_offset];
      flm[ind] +=
	ssign 
	* elfactor
	* expsm[m]
	* dl[0*dl_stride + m + dl_offset]
	* elssign * dl[0*dl_stride - spinneg + dl_offset]
	* Fmm[m*Fmm_stride + 0 + Fmm_offset];
//	* Gmm[(0+Fmm_offset)*Gmm_stride + m];
    }

    for (mm=1; mm<=el; mm++) {
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign 
	  * elfactor
	  * expsm[m]
	  * dl[mm*dl_stride + m + dl_offset]
	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
	  * (  Fmm[m*Fmm_stride + mm + Fmm_offset] 
	      + signs[m] * ssign
	       * Fmm[m*Fmm_stride - mm + Fmm_offset]);
	  /* * ( Gmm[(mm+Fmm_offset)*Gmm_stride + m] */
	  /*     + signs[m] * ssign */
	  /*     * Gmm[(-mm+Fmm_offset)*Gmm_stride + m]); */
      }

    }  

  }

  // Set flm values for negative m using conjugate symmetry.
  for (el=abs(spin); el<=L-1; el++) {
    for (m=1; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      ssht_sampling_elm2ind(&ind_nm, el, -m);
      flm[ind_nm] = signs[m] * conj(flm[ind]);
    }
  }

  // Free memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);
  free(Fmt);
  free(Fmm);
  free(inout);
  free(w);
  free(wr);
  free(Fmm_pad);
  free(tmp_pad);
  free(Gmm);
  free(sqrt_tbl);
  free(signs); 
  free(expsm);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s %s", SSHT_PROMPT, "Forward transform computed!");  

}
