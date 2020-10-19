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
#include <tgmath.h>  // Must be before fftw3.h
#include <fftw3.h>

#include "ssht/ssht_types.h"
#include "ssht/ssht_error.h"
#include "ssht/ssht_dl.h"
#include "ssht/ssht_sampling.h"
#include "ssht/ssht_adjoint.h"


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
void ssht_adjoint_mw_inverse_sov_sym(ssht_complex_double *flm, 
				     const ssht_complex_double *f, 
				     int L, int spin, 
				     ssht_dl_method_t dl_method,
				     int verbosity) {

  int el, m, mm, ind, t;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan;
  ssht_complex_double *inout;
  ssht_complex_double *Fmt, *Fmm;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *expsm, *expsmm;
  int exps_offset;
  int elmmsign, elssign;
  int spinneg;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  expsm = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsm)
  expsmm = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
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
    printf("%s%s\n", SSHT_PROMPT, 
	   "Computing adjoint inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT, 
	     "Using routine ssht_adjoint_mw_inverse_sov_sym...");
  }

  // Compute Fourier transform over phi, i.e. compute Fmt.
  Fmt = (ssht_complex_double*)calloc((2*L-1)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmt)
  Fmt_stride = 2*L-1;
  Fmt_offset = L-1;
  f_stride = 2*L-1;
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (t=0; t<=L-1; t++) {
    memcpy(inout, &f[t*f_stride], f_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(m=0; m<=L-1; m++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m];
    for(m=-(L-1); m<=-1; m++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m+2*L-1]; 
  }

  // Apply adjoint of periodic entension.
  for (m=-(L-1); m<=L-1; m++) 
    for (t=L; t<=2*L-2; t++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = 0.0;

  // Compute Fourier transform over theta, i.e. compute Fmm.
  Fmm = (ssht_complex_double*)calloc((2*L-1)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L-1;
  Fmm_offset = L-1;
  for (m=-(L-1); m<=L-1; m++) {
    memcpy(inout, &Fmt[(m+Fmt_offset)*Fmt_stride], Fmt_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L-1; mm++) 
      Fmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] = 
	inout[mm]; 
    for(mm=-(L-1); mm<=-1; mm++) 
      Fmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] = 
	inout[mm+2*L-1]; 
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (mm=-(L-1); mm<=L-1; mm++)
    for (m=-(L-1); m<=L-1; m++)
      Fmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] *=
  	expsmm[mm + exps_offset];

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
	* Fmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
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
	* Fmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
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
	  * ( Fmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]
      	      + signs[-m] * ssign
      	      * Fmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]);
      }
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign 
	  * elfactor
	  * expsm[m + exps_offset]
	  * dl[mm*dl_stride + m + dl_offset]
	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
	  * ( Fmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]
	      + signs[m] * ssign
	      * Fmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]);
      }

    }  

  }

  // Free memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);
  free(Fmt);
  free(Fmm);
  free(sqrt_tbl);
  free(signs);
  free(expsm);
  free(expsmm);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s%s", SSHT_PROMPT, "Adjoint inverse transform computed!");  

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
void ssht_adjoint_mw_inverse_sov_sym_real(ssht_complex_double *flm, 
					  const double *f, 
					  int L,
					  ssht_dl_method_t dl_method, 
					  int verbosity) {

  int el, m, mm, ind, ind_nm, t;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan;
  double *in_real;
  ssht_complex_double *inout, *out;
  ssht_complex_double *Fmt, *Fmm;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *expsm, *expsmm;
  int exps_offset;
  int elmmsign, elssign;
  int spinneg;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  expsm = (ssht_complex_double*)calloc(L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsm)
  expsmm = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
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
    printf("%s%s\n", SSHT_PROMPT, 
	   "Computing adjoint inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT, 
	     "Using routine ssht_adjoint_mw_inverse_sov_sym_real...");
  }

  // Compute Fourier transform over phi, i.e. compute Fmt.
  Fmt = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmt)
  Fmt_stride = 2*L-1;  
  f_stride = 2*L-1;
  in_real = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in_real)
  out = (ssht_complex_double*)calloc(L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out)
  plan = fftw_plan_dft_r2c_1d(2*L-1, in_real, out, FFTW_MEASURE);
  for (t=0; t<=L-1; t++) {
    memcpy(in_real, &f[t*f_stride], f_stride*sizeof(double));
        fftw_execute_dft_r2c(plan, in_real, out);
    for(m=0; m<=L-1; m++) 
      Fmt[m*Fmt_stride + t] = out[m];
  }
  free(in_real);
  free(out);
  fftw_destroy_plan(plan);

  // Apply adjoint of periodic extension.
  for (m=0; m<=L-1; m++) 
    for (t=L; t<=2*L-2; t++) 
      Fmt[m*Fmt_stride + t] = 0.0;

  // Compute Fourier transform over theta, i.e. compute Fmm.
  Fmm = (ssht_complex_double*)calloc((2*L-1)*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = L;
  Fmm_offset = L-1;
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (m=0; m<=L-1; m++) {
    memcpy(inout, &Fmt[m*Fmt_stride], Fmt_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L-1; mm++) 
      Fmm[(mm+Fmm_offset)*Fmm_stride + m] = 
	inout[mm]; 
    for(mm=-(L-1); mm<=-1; mm++) 
      Fmm[(mm+Fmm_offset)*Fmm_stride + m] = 
	inout[mm+2*L-1]; 
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (mm=-(L-1); mm<=L-1; mm++)
    for (m=0; m<=L-1; m++)
      Fmm[(mm+Fmm_offset)*Fmm_stride + m] *= 
	expsmm[mm + exps_offset];

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
	* Fmm[(0+Fmm_offset)*Fmm_stride + m];
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
  	  * ( Fmm[(mm+Fmm_offset)*Fmm_stride + m]
	      + signs[m] * ssign
	      * Fmm[(-mm+Fmm_offset)*Fmm_stride + m]);
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
  free(sqrt_tbl);
  free(signs); 
  free(expsm);
  free(expsmm);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s%s", SSHT_PROMPT, "Adjoint inverse transform computed!");  

}


/*!  
 * Compute adjoint of forward transform for MW method using separation
 * of variables, fast Fourier transforms, performing convolution with
 * weights as product in transformed space and exploiting all
 * symmetries (for complex spin signal).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_forward_sov_sym(ssht_complex_double *f, const ssht_complex_double *flm,
				     int L, int spin,
				     ssht_dl_method_t dl_method,
				     int verbosity) {

  int el, m, mm, ind;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  ssht_complex_double mmfactor;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *exps;
  int exps_offset;
  double elmmsign, elssign;
  int spinneg;
  ssht_complex_double *Fmm;
  int Fmm_offset, Fmm_stride; 

  fftw_plan plan, plan_bwd, plan_fwd;
  ssht_complex_double *Ftm, *Gmm;
  ssht_complex_double *w, *wr;
  int w_offset;
  ssht_complex_double *Fmm_pad, *tmp_pad;
  int f_stride, Ftm_stride, Ftm_offset;
  int r, t, p;
  ssht_complex_double *inout;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  exps = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(exps)
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
    exps[m + exps_offset] = cexp(-I*SSHT_PION2*(m+spin));

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing adjoint forward transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_adjoint_mw_forward_sov_sym...");
  }

  // Compute Fmm.
  Fmm = (ssht_complex_double*)calloc((2*L-1)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_offset = L-1;
  Fmm_stride = 2*L-1;
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
  inds_offset = L-1;
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

    // Compute Fmm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el *el + el;
    for (m=-el; m<=el; m++)
      inds[m + inds_offset] = el2pel + m;
    for (mm=0; mm<=el; mm++) {
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      for (m=-el; m<=-1; m++) {
	ind = inds[m + inds_offset];
    	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] +=
    	  ssign
    	  * elfactor
	  * exps[m + exps_offset]
    	  * elmmsign * dl[mm*dl_stride - m + dl_offset]
    	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
    	  * flm[ind];
      }
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
    	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] +=
    	  ssign
    	  * elfactor
	  * exps[m + exps_offset]
    	  * dl[mm*dl_stride + m + dl_offset]
    	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
    	  * flm[ind];
      }

    }

  }

  // Free dl memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);

  // Use symmetry to compute Fmm for negative mm.
  for (m=-(L-1); m<=L-1; m++)
    for (mm=-(L-1); mm<=-1; mm++)
      Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] =
	signs[abs(m)] * ssign
	* Fmm[(m + Fmm_offset)*Fmm_stride - mm + Fmm_offset];

  // Compute weights.
  w = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(w)
  w_offset = 2*(L-1);
  for (mm=-2*(L-1); mm<=2*(L-1); mm++)
    w[mm+w_offset] = ssht_sampling_weight_mw(mm);

  // Compute IFFT of w to give wr.
  wr = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(wr)
  inout = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
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
  Fmm_pad = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_pad)
  tmp_pad = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(tmp_pad)
  Gmm = (ssht_complex_double*)calloc((2*L-1)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Gmm)
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
      Gmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] =
	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);

  }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (mm=-(L-1); mm<=L-1; mm++) {
    mmfactor = cexp(I*mm*SSHT_PI/(2.0*L-1.0));
    for (m=-(L-1); m<=L-1; m++)
      Gmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] *=
  	mmfactor;
  }

  // Compute Fourier transform over theta.
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  Ftm = (ssht_complex_double*)calloc((2*L-1)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = 2*L-1;
  Ftm_offset = L-1;
  for (m=-(L-1); m<=L-1; m++) {

    for(mm=0; mm<=L-1; mm++)
      inout[mm] = Gmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset];
    for(mm=-(L-1); mm<=-1; mm++)
      inout[mm+2*L-1] = Gmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset];
    fftw_execute_dft(plan, inout, inout);

    for(t=0; t<=2*L-2; t++)
      Ftm[t*Ftm_stride + m + Ftm_offset] = inout[t] / (2.0*L-1.0);

  }

  // Adjoint of periodic extension of Ftm.
  for(t=0; t<=L-2; t++)
    for (m=-(L-1); m<=L-1; m++)
      Ftm[t*Ftm_stride + m + Ftm_offset] = 
	Ftm[t*Ftm_stride + m + Ftm_offset] 
	+ signs[abs(m)] * ssign * Ftm[(2*L-2-t)*Ftm_stride + m + Ftm_offset];

  // Compute Fourier transform over phi.
  f_stride = 2*L-1;
  for(t=0; t<=L-1; t++) {

    for(m=0; m<=L-1; m++)
      inout[m] = Ftm[t*Ftm_stride + m + Ftm_offset];
    for(m=-(L-1); m<=-1; m++)
      inout[m+2*L-1] = Ftm[t*Ftm_stride + m + Ftm_offset];
    fftw_execute_dft(plan, inout, inout);

    for(p=0; p<=2*L-2; p++)
      f[t*f_stride + p] = inout[p] / (2.0*L-1.0);

  }
  fftw_destroy_plan(plan);
  free(inout);

  // Free memory.
  free(Fmm);
  free(Ftm);
  free(w);
  free(wr);
  free(Fmm_pad);
  free(tmp_pad);
  free(Gmm);

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Adjoint forward transform computed!");

}


/*!  
 * Compute adjoint of forward transform for MW method of real scalar
 * signal using separation of variables, fast Fourier transforms,
 * performing convolution with weights as product in transformed space
 * and exploiting all symmetries (including additional symmetries for
 * real signals).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_forward_sov_sym_real(double *f, 
					  const ssht_complex_double *flm,
					  int L,
					  ssht_dl_method_t dl_method,
					  int verbosity) {

  int el, m, mm, ind;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  ssht_complex_double mmfactor;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *exps;
  int exps_offset;
  double elmmsign, elssign;
  int spinneg;
  ssht_complex_double *Fmm;
  int Fmm_offset, Fmm_stride;

  fftw_plan plan, plan_bwd, plan_fwd;
  ssht_complex_double *Ftm, *Gmm;
  ssht_complex_double *w, *wr;
  int w_offset;
  ssht_complex_double *Fmm_pad, *tmp_pad;
  int f_stride, Ftm_stride;
  int r, t, p;
  ssht_complex_double *inout;
  ssht_complex_double *in;
  double *out_real;
  int Gmm_stride;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  exps = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(exps)
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
    exps[m + exps_offset] = cexp(-I*SSHT_PION2*(m+spin));

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing adjoint forward transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_adjoint_mw_forward_sov_sym_real...");
  }

  // Compute Fmm.
  Fmm = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_offset = L-1;
  Fmm_stride = 2*L-1;    
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);   
  inds_offset = L-1;
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

    // Compute Fmm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el *el + el;    
    for (m=0; m<=el; m++)
      inds[m + inds_offset] = el2pel + m; 
    for (mm=0; mm<=el; mm++) {
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
    	Fmm[m*Fmm_stride + mm + Fmm_offset] +=
    	  ssign
    	  * elfactor
	  * exps[m + exps_offset]    	  
    	  * dl[mm*dl_stride + m + dl_offset]
    	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
    	  * flm[ind];
      }

    }

  }

  // Free dl memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);

  // Use symmetry to compute Fmm for negative mm.
  for (m=0; m<=L-1; m++) 
    for (mm=-(L-1); mm<=-1; mm++) 
      Fmm[m*Fmm_stride + mm + Fmm_offset] = 
	signs[abs(m)] * ssign 
	* Fmm[m*Fmm_stride - mm + Fmm_offset];

  // Compute weights.
  w = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(w)
  w_offset = 2*(L-1);
  for (mm=-2*(L-1); mm<=2*(L-1); mm++)
    w[mm+w_offset] = ssht_sampling_weight_mw(mm);

  // Compute IFFT of w to give wr.
  wr = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(wr)
  inout = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
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
  Fmm_pad = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_pad)
  tmp_pad = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(tmp_pad)
  Gmm = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  Gmm_stride = 2*L-1;
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
      Gmm[m*Gmm_stride + mm + Fmm_offset] = 
	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);

  }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (mm=-(L-1); mm<=L-1; mm++) {
    mmfactor = cexp(I*mm*SSHT_PI/(2.0*L-1.0));
    for (m=0; m<=L-1; m++) 
      Gmm[m*Gmm_stride + mm + Fmm_offset] *= 
	mmfactor;
  }

  // Compute Fourier transform over theta.
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  Ftm = (ssht_complex_double*)calloc((2*L-1)*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = L;
  for (m=0; m<=L-1; m++) {

    for(mm=0; mm<=L-1; mm++)
      inout[mm] = Gmm[m*Gmm_stride + mm + Fmm_offset];
    for(mm=-(L-1); mm<=-1; mm++)
      inout[mm+2*L-1] = Gmm[m*Gmm_stride + mm + Fmm_offset];
    fftw_execute_dft(plan, inout, inout);

    for(t=0; t<=2*L-2; t++)
      Ftm[t*Ftm_stride + m] = inout[t] / (2.0*L-1.0);

  }
  fftw_destroy_plan(plan);
  free(inout);

  // Adjoint of periodic extension of Ftm.
  for(t=0; t<=L-2; t++)
    for (m=0; m<=L-1; m++)
      Ftm[t*Ftm_stride + m] =
  	Ftm[t*Ftm_stride + m]
  	+ signs[abs(m)] * ssign * Ftm[(2*L-2-t)*Ftm_stride + m];

  // Compute Fourier transform over phi.
  out_real = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out_real)
  in = (ssht_complex_double*)calloc(L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in)
  plan = fftw_plan_dft_c2r_1d(2*L-1, in, out_real, FFTW_MEASURE);
  f_stride = 2*L-1;
  for(t=0; t<=L-1; t++) {
    
    memcpy(in, &Ftm[t*Ftm_stride], Ftm_stride*sizeof(ssht_complex_double));
    fftw_execute_dft_c2r(plan, in, out_real);

    for(p=0; p<=2*L-2; p++)
      f[t*f_stride + p] = out_real[p] / (2.0*L-1.0);

  }
  fftw_destroy_plan(plan);
  free(out_real);
  free(in);

  // Free memory.
  free(Fmm);
  free(Ftm);
  free(w);
  free(wr);
  free(Fmm_pad);
  free(tmp_pad);
  free(Gmm);

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Adjoint forward transform computed!");

}


//============================================================================
// MW South pole interfaces
//============================================================================


/*!  
 * South pole wrapper for adjoint forward transform for MW method.
 * The South pole is defined by a single sample and its corresponding
 * phi angle, rather than specifying samples for all phi at the South
 * pole (which are simply related by the rotation of a spin function
 * in its tangent plane).
 *
 * \param[out] f Function on sphere (excluding South pole).
 * \param[out] f_sp Function sample on South pole.
 * \param[out] phi_sp Phi angle corresponding to quoted sample at 
 * South pole.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_forward_sov_sym_pole(ssht_complex_double *f, 
					  ssht_complex_double *f_sp, double *phi_sp,
					  ssht_complex_double *flm, 
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity) {

  ssht_complex_double* f_full;
  int f_stride = 2*L-1;

  // Allocate full array.
  f_full = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)

  // Perform inverse transform.
  ssht_adjoint_mw_forward_sov_sym(f_full, flm, L, spin, 
				  dl_method, verbosity);	  

  // Copy output function values, including separate point for South pole.
  memcpy(f, f_full, (L-1)*(2*L-1)*sizeof(ssht_complex_double));
  *f_sp = f_full[(L-1)*f_stride + 0];
  *phi_sp = ssht_sampling_mw_p2phi(0, L);
	
  // Free memory.
  free(f_full);

}


/*!  
 * South pole wrapper for adjoint forward transform of real scalar
 * function for MW method.  The South pole is defined by a single
 * sample, rather than specifying samples for all phi at the South
 * pole (which for a scalar function are identical).
 *
 * \param[out] f Function on sphere (excluding South pole).
 * \param[out] f_sp Function sample on South pole.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_forward_sov_sym_real_pole(double *f, 
					       double *f_sp,
					       ssht_complex_double *flm, 
					       int L, 
					       ssht_dl_method_t dl_method, 
					       int verbosity) {

  double *f_full;
  int f_stride = 2*L-1;

  // Allocate full array.
  f_full = (double*)calloc(L*(2*L-1), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)

  // Perform inverse transform.
  ssht_adjoint_mw_forward_sov_sym_real(f_full, flm, L, 
				       dl_method, verbosity);

  // Copy output function values, including separate point for South pole.
  memcpy(f, f_full, (L-1)*(2*L-1)*sizeof(double));
  *f_sp = f_full[(L-1)*f_stride + 0];

  // Free memory.
  free(f_full);

}


/*!  
 * South pole wrapper for adjoint inverse transform for MW method.
 * The South pole is defined by a single sample and its corresponding
 * phi angle, rather than specifying samples for all phi at the South
 * pole (which are simply related by the rotation of a spin function
 * in its tangent plane).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere (excluding South pole).
 * \param[in] f_sp Function sample on South pole.
 * \param[in] phi_sp Phi angle corresponding to quoted sample at 
 * South pole.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_inverse_sov_sym_pole(ssht_complex_double *flm, ssht_complex_double *f,
					  ssht_complex_double f_sp, double phi_sp,
					  int L, int spin, 
					  ssht_dl_method_t dl_method,
					  int verbosity) {

  ssht_complex_double *f_full;
  int p, f_stride = 2*L-1;
  double phi;

  // Copy function values to full array.
  f_full = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)
  memcpy(f_full, f, (L-1)*(2*L-1)*sizeof(ssht_complex_double));

  // Define South pole for all phi.
  for (p=0; p<=2*L-2; p++) {
    phi = ssht_sampling_mw_p2phi(p, L);
    f_full[(L-1)*f_stride + p] = f_sp * cexp(I*spin*(phi-phi_sp)); 
  }

  // Perform forward transform.
  ssht_adjoint_mw_inverse_sov_sym(flm, f_full, L, spin, 
				  dl_method, verbosity);

  // Free memory.
  free(f_full);

}


/*!  
 * South pole wrapper for adjoint inverse transform of real scalar
 * function for MW method.  The South pole is defined by a single
 * sample, rather than specifying samples for all phi at the South
 * pole (which for a scalar function are identical).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere (excluding South pole).
 * \param[in] f_sp Function sample on South pole.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_inverse_sov_sym_real_pole(ssht_complex_double *flm, 
					       double *f, 
					       double f_sp,
					       int L, 
					       ssht_dl_method_t dl_method,
					       int verbosity) {

  double *f_full;
  int p, f_stride = 2*L-1;

  // Copy function values to full array.
  f_full = (double*)calloc(L*(2*L-1), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)
  memcpy(f_full, f, (L-1)*(2*L-1)*sizeof(double));

  // Define South pole for all phi.
  for (p=0; p<=2*L-2; p++)
    f_full[(L-1)*f_stride + p] = f_sp; 

  // Perform forward transform.
  ssht_adjoint_mw_inverse_sov_sym_real(flm, f_full, L, 
				       dl_method, verbosity);

  // Free memory.
  free(f_full);

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
void ssht_adjoint_mw_inverse_sov_sym_ss(ssht_complex_double *flm, ssht_complex_double *f, 
					int L, int spin, 
					ssht_dl_method_t dl_method,
					int verbosity) {

  int el, m, mm, ind, t;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan;
  ssht_complex_double *inout;
  ssht_complex_double *Fmt, *Fmm;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *expsm, *expsmm;
  int exps_offset;
  int elmmsign, elssign;
  int spinneg;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  expsm = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(expsm)
  expsmm = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
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
    printf("%s%s\n", SSHT_PROMPT, 
	   "Computing adjoint inverse transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT, 
	     "Using routine ssht_adjoint_mw_inverse_sov_sym_ss...");
  }

  // Compute Fourier transform over phi, i.e. compute Fmt.
  // Note that t and p indices of fext are increased in size by
  // one compared to usual sampling.
  Fmt = (ssht_complex_double*)calloc((2*L)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmt)
  Fmt_stride = 2*L;
  Fmt_offset = L-1;
  f_stride = 2*L;
  inout = (ssht_complex_double*)calloc(2*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (t=0; t<=L; t++) {
    memcpy(inout, &f[t*f_stride], f_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(m=0; m<=L; m++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m]; 
    for(m=-(L-1); m<=-1; m++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m+2*L-1+1]; 
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Apply adjoint of periodic extension.
  for (m=-(L-1); m<=L; m++) 
    for (t=L+1; t<=2*L-1; t++) 
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = 0.0;

  // Compute Fourier transform over theta, i.e. compute Fmm.
  // Note that m and mm indices are increased in size by one.
  Fmm = (ssht_complex_double*)calloc((2*L)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L;
  Fmm_offset = L-1;
  inout = (ssht_complex_double*)calloc(2*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (m=-(L-1); m<=L; m++) {
    memcpy(inout, &Fmt[(m+Fmt_offset)*Fmt_stride], Fmt_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L; mm++) 
      Fmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] = 
	inout[mm]; 
    for(mm=-(L-1); mm<=-1; mm++) 
      Fmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] = 
	inout[mm+2*L-1+1]; 
  }
  fftw_destroy_plan(plan);
  free(inout);

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
	* Fmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
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
	* Fmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
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
      	  * ( Fmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]
      	      + signs[-m] * ssign
      	      * Fmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]);
      }
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign 
	  * elfactor
	  * expsm[m + exps_offset]
	  * dl[mm*dl_stride + m + dl_offset]
	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
	  * ( Fmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]
	      + signs[m] * ssign
	      * Fmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]);
      }

    }  

  }

  // Free memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);
  free(Fmt);
  free(Fmm);
  free(sqrt_tbl);
  free(signs); 
  free(expsm);
  free(expsmm);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s%s", SSHT_PROMPT, "Adjoint inverse transform computed!");  

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
void ssht_adjoint_mw_inverse_sov_sym_ss_real(ssht_complex_double *flm, double *f, 
					     int L, 
					     ssht_dl_method_t dl_method, 
					     int verbosity) {

  int el, m, mm, ind, ind_nm, t;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan;
  double *in_real;
  ssht_complex_double *inout, *out;
  ssht_complex_double *Fmt, *Fmm;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *expsm;
  int exps_offset;
  int elmmsign, elssign;
  int spinneg;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  expsm = (ssht_complex_double*)calloc(L, sizeof(ssht_complex_double));
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
  for (m=0; m<=L-1; m++)
    expsm[m] = cexp(I*SSHT_PION2*(m+spin));

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT, 
	   "Computing adjoint inverse transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT, 
	     "Using routine ssht_adjoint_mw_inverse_sov_sym_ss_real...");
  }

  // Compute Fourier transform over phi, i.e. compute Fmt.
  // Note that t and p indices of fext are increased in size by
  // one compared to usual sampling.
  Fmt = (ssht_complex_double*)calloc((L+1)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmt)
  Fmt_stride = 2*L;
  f_stride = 2*L;
  in_real = (double*)calloc(2*L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in_real)
  out = (ssht_complex_double*)calloc(L+1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out)
  plan = fftw_plan_dft_r2c_1d(2*L, in_real, out, FFTW_MEASURE);
  for (t=0; t<=L; t++) {
    memcpy(in_real, &f[t*f_stride], f_stride*sizeof(double));
        fftw_execute_dft_r2c(plan, in_real, out);
    for(m=0; m<=L; m++) 
      Fmt[m*Fmt_stride + t] = out[m]; 
  }
  free(in_real);
  free(out);
  fftw_destroy_plan(plan);

  // Apply adjoint of periodic extension.
  for (m=0; m<=L; m++) 
    for (t=L+1; t<=2*L-1; t++) 
      Fmt[m*Fmt_stride + t] = 0.0;

  // Compute Fourier transform over theta, i.e. compute Fmm.
  // Note that m and mm indices are increased in size by one.
  Fmm = (ssht_complex_double*)calloc((2*L)*(L+1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = L+1;
  Fmm_offset = L-1;
  inout = (ssht_complex_double*)calloc(2*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (m=0; m<=L; m++) {
    memcpy(inout, &Fmt[m*Fmt_stride], Fmt_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L; mm++) 
      Fmm[(mm+Fmm_offset)*Fmm_stride + m] = 
	inout[mm]; 
    for(mm=-(L-1); mm<=-1; mm++) 
      Fmm[(mm+Fmm_offset)*Fmm_stride + m] =
	inout[mm+2*L-1+1]; 
  }
  fftw_destroy_plan(plan);
  free(inout);

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
	* Fmm[(0+Fmm_offset)*Fmm_stride + m];
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
	  * ( Fmm[(mm+Fmm_offset)*Fmm_stride + m]
	      + signs[m] * ssign
	      * Fmm[(-mm+Fmm_offset)*Fmm_stride + m]);
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
  free(sqrt_tbl);
  free(signs); 
  free(expsm);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s%s", SSHT_PROMPT, "Adjoint inverse transform computed!");  

}


/*!  
 * Compute adjoint of forward transform for MW method with symmetric
 * sampling using separation of variables, fast Fourier transforms,
 * performing convolution with weights as product in transformed space
 * and exploiting all symmetries (for complex spin signal).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_forward_sov_sym_ss(ssht_complex_double *f, ssht_complex_double *flm,
					int L, int spin,
					ssht_dl_method_t dl_method,
					int verbosity) {

  int el, m, mm, ind;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *exps;
  int exps_offset;
  double elmmsign, elssign;
  int spinneg;
  ssht_complex_double *Fmm;
  int Fmm_offset, Fmm_stride;

  fftw_plan plan, plan_bwd, plan_fwd;
  ssht_complex_double *Ftm, *Gmm;
  ssht_complex_double *w, *wr;
  int w_offset;
  ssht_complex_double *Fmm_pad, *tmp_pad;
  int f_stride, Ftm_stride, Ftm_offset;
  int r, t, p;
  ssht_complex_double *inout;

  int Gmm_offset, Gmm_stride;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  exps = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(exps)
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
    exps[m + exps_offset] = cexp(-I*SSHT_PION2*(m+spin));

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing adjoint forward transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_adjoint_mw_forward_sov_sym_ss...");
  }

  // Compute Fmm.
  Fmm = (ssht_complex_double*)calloc((2*L-1)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_offset = L-1;
  Fmm_stride = 2*L-1;
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
  inds_offset = L-1;
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

    // Compute Fmm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el *el + el;
    for (m=-el; m<=el; m++)
      inds[m + inds_offset] = el2pel + m;
    for (mm=0; mm<=el; mm++) {
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      for (m=-el; m<=-1; m++) {
	ind = inds[m + inds_offset];
    	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] +=
    	  ssign
    	  * elfactor
	  * exps[m + exps_offset]
    	  * elmmsign * dl[mm*dl_stride - m + dl_offset]
    	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
    	  * flm[ind];
      }
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
    	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] +=
    	  ssign
    	  * elfactor
	  * exps[m + exps_offset]
    	  * dl[mm*dl_stride + m + dl_offset]
    	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
    	  * flm[ind];
      }

    }

  }

  // Free dl memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);

  // Use symmetry to compute Fmm for negative mm.
  for (m=-(L-1); m<=L-1; m++)
    for (mm=-(L-1); mm<=-1; mm++)
      Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] =
	signs[abs(m)] * ssign
	* Fmm[(m + Fmm_offset)*Fmm_stride - mm + Fmm_offset];

  // Compute weights.
  w = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(w)
  w_offset = 2*(L-1);
  for (mm=-2*(L-1); mm<=2*(L-1); mm++)
    w[mm+w_offset] = ssht_sampling_weight_mw(mm);

  // Compute IFFT of w to give wr.
  wr = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(wr)
  inout = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
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
  Fmm_pad = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_pad)
  tmp_pad = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(tmp_pad)
  Gmm = (ssht_complex_double*)calloc((2*L)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Gmm)
  Gmm_stride = 2*L;
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
      Gmm[(m+Gmm_offset)*Gmm_stride + mm + Gmm_offset] =
	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);

  }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);
  free(inout);

  // Compute Fourier transform over theta.
  inout = (ssht_complex_double*)calloc(2*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  Ftm = (ssht_complex_double*)calloc((2*L)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = 2*L;
  Ftm_offset = L-1;
  for (m=-(L-1); m<=L; m++) {

    for(mm=0; mm<=L; mm++)
      inout[mm] = Gmm[(m+Gmm_offset)*Gmm_stride + mm + Gmm_offset];
    for(mm=-(L-1); mm<=-1; mm++)
      inout[mm+2*L-1+1] = Gmm[(m+Gmm_offset)*Gmm_stride + mm + Gmm_offset];
    fftw_execute_dft(plan, inout, inout);

    for(t=0; t<=2*L-1; t++)
      Ftm[t*Ftm_stride + m + Ftm_offset] = inout[t] / (2.0*L);

  }

  // Adjoint of periodic extension of Ftm.
  for(t=1; t<=L-1; t++)
    for (m=-(L-1); m<=L; m++)
      Ftm[t*Ftm_stride + m + Ftm_offset] = 
	Ftm[t*Ftm_stride + m + Ftm_offset] 
	+ signs[abs(m)] * ssign * Ftm[(2*L-t)*Ftm_stride + m + Ftm_offset];

  // Compute Fourier transform over phi.
  f_stride = 2*L;
  for(t=0; t<=L; t++) {

    for(m=0; m<=L; m++)
      inout[m] = Ftm[t*Ftm_stride + m + Ftm_offset];
    for(m=-(L-1); m<=-1; m++)
      inout[m+2*L-1+1] = Ftm[t*Ftm_stride + m + Ftm_offset];
    fftw_execute_dft(plan, inout, inout);

    for(p=0; p<=2*L-1; p++)
      f[t*f_stride + p] = inout[p] / (2.0*L);

  }
  fftw_destroy_plan(plan);
  free(inout);

  // Free memory.
  free(Fmm);
  free(Ftm);
  free(w);
  free(wr);
  free(Fmm_pad);
  free(tmp_pad);
  free(Gmm);

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Adjoint forward transform computed!");

}


/*!  
 * Compute adjoint of forward transform for MW method with symmetric
 * sampling of real scalar signal using separation of variables, fast
 * Fourier transforms, performing convolution with weights as product
 * in transformed space and exploiting all symmetries (including
 * additional symmetries for real signals).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_forward_sov_sym_ss_real(double *f, 
					     ssht_complex_double *flm,
					     int L,
					     ssht_dl_method_t dl_method,
					     int verbosity) {

  int el, m, mm, ind;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *exps;
  int exps_offset;
  double elmmsign, elssign;
  int spinneg;
  ssht_complex_double *Fmm;
  int Fmm_offset, Fmm_stride;

  fftw_plan plan, plan_bwd, plan_fwd;
  ssht_complex_double *Ftm, *Gmm;
  ssht_complex_double *w, *wr;
  int w_offset;
  ssht_complex_double *Fmm_pad, *tmp_pad;
  int f_stride, Ftm_stride;
  int r, t, p;
  ssht_complex_double *inout;
  ssht_complex_double *in;
  double *out_real;
  int Gmm_offset, Gmm_stride;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  exps = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(exps)
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
    exps[m + exps_offset] = cexp(-I*SSHT_PION2*(m+spin));

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing adjoint forward transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_adjoint_mw_forward_sov_sym_ss_real...");
  }

  // Compute Fmm.
  Fmm = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_offset = L-1;
  Fmm_stride = 2*L-1;    
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);   
  inds_offset = L-1;
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

    // Compute Fmm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el *el + el;    
    for (m=0; m<=el; m++)
      inds[m + inds_offset] = el2pel + m; 
    for (mm=0; mm<=el; mm++) {
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
    	Fmm[m*Fmm_stride + mm + Fmm_offset] +=
    	  ssign
    	  * elfactor
	  * exps[m + exps_offset]    	  
    	  * dl[mm*dl_stride + m + dl_offset]
    	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
    	  * flm[ind];
      }

    }

  }

  // Free dl memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);

  // Use symmetry to compute Fmm for negative mm.
  for (m=0; m<=L-1; m++) 
    for (mm=-(L-1); mm<=-1; mm++) 
      Fmm[m*Fmm_stride + mm + Fmm_offset] =      
	signs[abs(m)] * ssign 
	* Fmm[m*Fmm_stride - mm + Fmm_offset];

  // Compute weights.
  w = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(w)
  w_offset = 2*(L-1);
  for (mm=-2*(L-1); mm<=2*(L-1); mm++)
    w[mm+w_offset] = ssht_sampling_weight_mw(mm);

  // Compute IFFT of w to give wr.
  wr = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(wr)
  inout = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
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
  Fmm_pad = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_pad)
  tmp_pad = (ssht_complex_double*)calloc(4*L-3, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(tmp_pad)
  Gmm = (ssht_complex_double*)calloc((L+1)*(2*L), sizeof(ssht_complex_double));
  Gmm_stride = 2*L;
  Gmm_offset = L-1;
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
      Gmm[m*Gmm_stride + mm + Gmm_offset] = 
	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);

  }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);
  free(inout);
  
  // Compute Fourier transform over theta.
  inout = (ssht_complex_double*)calloc(2*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  Ftm = (ssht_complex_double*)calloc((2*L)*(L+1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = L+1;
  for (m=0; m<=L; m++) {

    for(mm=0; mm<=L; mm++)
      inout[mm] = Gmm[m*Gmm_stride + mm + Gmm_offset];
    for(mm=-(L-1); mm<=-1; mm++)
      inout[mm+2*L-1+1] = Gmm[m*Gmm_stride + mm + Gmm_offset];
    fftw_execute_dft(plan, inout, inout);

    for(t=0; t<=2*L-1; t++)
      Ftm[t*Ftm_stride + m] = inout[t] / (2.0*L);

  }
  fftw_destroy_plan(plan);
  free(inout);

  // Adjoint of periodic extension of Ftm.
  for(t=1; t<=L-1; t++)
    for (m=0; m<=L; m++)
      Ftm[t*Ftm_stride + m] =
  	Ftm[t*Ftm_stride + m]
  	+ signs[abs(m)] * ssign * Ftm[(2*L-t)*Ftm_stride + m];

  // Compute Fourier transform over phi.
  out_real = (double*)calloc(2*L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out_real)
  in = (ssht_complex_double*)calloc(L+1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in)
  plan = fftw_plan_dft_c2r_1d(2*L, in, out_real, FFTW_MEASURE);
  f_stride = 2*L;
  for(t=0; t<=L; t++) {
    
    memcpy(in, &Ftm[t*Ftm_stride], Ftm_stride*sizeof(ssht_complex_double));
    fftw_execute_dft_c2r(plan, in, out_real);

    for(p=0; p<=2*L-1; p++)
      f[t*f_stride + p] = out_real[p] / (2.0*L);

  }
  fftw_destroy_plan(plan);
  free(out_real);
  free(in);

  // Free memory.
  free(Fmm);
  free(Ftm);
  free(w);
  free(wr);
  free(Fmm_pad);
  free(tmp_pad);
  free(Gmm);

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Adjoint forward transform computed!");

}


//============================================================================
// MW SS Noth-South pole interfaces
//============================================================================


/*!  
 * North-South pole wrapper for adjoint forward transform for MW
 * method with symmetric sampling.  The poles are defined by single
 * samples and their corresponding phi angle, rather than specifying
 * samples for all phi at the poles (which are simply related by the
 * rotation of a spin function in its tangent plane).
 *
 * \param[out] f Function on sphere (excluding poles).
 * \param[out] f_np Function sample on North pole.
 * \param[out] phi_np Phi angle corresponding to quoted sample at 
 * North pole.
 * \param[out] f_sp Function sample on South pole.
 * \param[out] phi_sp Phi angle corresponding to quoted sample at 
 * South pole.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_forward_sov_sym_ss_pole(ssht_complex_double *f, 
					     ssht_complex_double *f_np, double *phi_np,
					     ssht_complex_double *f_sp, double *phi_sp,
					     ssht_complex_double *flm, 
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity) {

  ssht_complex_double* f_full;
  int t, f_stride = 2*L;

  // Allocate full array.
  f_full = (ssht_complex_double*)calloc((L+1)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)

  // Perform inverse transform.
  ssht_adjoint_mw_forward_sov_sym_ss(f_full, flm, L, spin, 
				     dl_method, verbosity);	  

  // Copy output function values, including separate points for  poles.
 for (t=1; t<=L-1; t++)
   memcpy(&f[(t-1)*f_stride], &f_full[t*f_stride], 
	  (2*L)*sizeof(ssht_complex_double));
  *f_np = f_full[0];
  *phi_np = ssht_sampling_mw_ss_p2phi(0, L);
  *f_sp = f_full[L*f_stride + 0];
  *phi_sp = ssht_sampling_mw_ss_p2phi(0, L);
	
  // Free memory.
  free(f_full);

}


/*!  
 * North-South pole wrapper for adjoint forward transform of real
 * scalar function for MW method with symmetric sampling.  The poles
 * are defined by single samples, rather than specifying samples for
 * all phi at the poles (which for a scalar function are identical).
 *
 * \param[out] f Function on sphere (excluding poles).
 * \param[out] f_sp Function sample on South pole.
 * \param[out] f_np Function sample on North pole.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_forward_sov_sym_ss_real_pole(double *f, 
						  double *f_np,
						  double *f_sp,
						  ssht_complex_double *flm, 
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity) {

  double *f_full;
  int t, f_stride = 2*L;

  // Allocate full array.
  f_full = (double*)calloc((L+1)*(2*L), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)

  // Perform inverse transform.
  ssht_adjoint_mw_forward_sov_sym_ss_real(f_full, flm, L, 
					  dl_method, verbosity);

  // Copy output function values, including separate points for  poles.
  for (t=1; t<=L-1; t++)
   memcpy(&f[(t-1)*f_stride], &f_full[t*f_stride], 
	  (2*L)*sizeof(double));
  *f_np = f_full[0];
  *f_sp = f_full[L*f_stride + 0];

  // Free memory.
  free(f_full);

}


/*!  
 * North-South pole wrapper for adjoint inverse transform for MW
 * method with symmetric sampling.  The poles are defined by single
 * samples and their corresponding phi angle, rather than specifying
 * samples for all phi at the poles (which are simply related by the
 * rotation of a spin function in its tangent plane).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere (excluding poles).
 * \param[in] f_np Function sample on North pole.
 * \param[in] phi_np Phi angle corresponding to quoted sample at 
 * North pole.
 * \param[in] f_sp Function sample on South pole.
 * \param[in] phi_sp Phi angle corresponding to quoted sample at 
 * South pole.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_inverse_sov_sym_ss_pole(ssht_complex_double *flm, ssht_complex_double *f,
					     ssht_complex_double f_np, double phi_np,
					     ssht_complex_double f_sp, double phi_sp,
					     int L, int spin, 
					     ssht_dl_method_t dl_method,
					     int verbosity) {

  ssht_complex_double *f_full;
  int t, p, f_stride = 2*L;
  double phi;

  // Copy function values to full array.
  f_full = (ssht_complex_double*)calloc((L+1)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)
  for (t=1; t<=L-1; t++)
    memcpy(&f_full[t*f_stride], &f[(t-1)*f_stride], 
	   (2*L)*sizeof(ssht_complex_double));

  // Define poles for all phi.
  for (p=0; p<=2*L-1; p++) {
    phi = ssht_sampling_mw_ss_p2phi(p, L);
    f_full[0*f_stride + p] = f_np * cexp(-I*spin*(phi-phi_np)); 
    f_full[L*f_stride + p] = f_sp * cexp(I*spin*(phi-phi_sp)); 
  }

  // Perform forward transform.
  ssht_adjoint_mw_inverse_sov_sym_ss(flm, f_full, L, spin, 
				     dl_method, verbosity);

  // Free memory.
  free(f_full);

}


/*!  
 * North-South pole wrapper for adjoint inverse transform of real
 * scalar function for MW method with symmetric sampling.  The poles
 * are defined by single samples, rather than specifying samples for
 * all phi at the poles (which for a scalar function are identical).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere (excluding poles).
 * \param[in] f_np Function sample on North pole.
 * \param[in] f_sp Function sample on South pole.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_adjoint_mw_inverse_sov_sym_ss_real_pole(ssht_complex_double *flm, 
						  double *f, 
						  double f_np,
						  double f_sp,
						  int L, 
						  ssht_dl_method_t dl_method,
						  int verbosity) {

  double *f_full;
  int t, p, f_stride = 2*L;

  // Copy function values to full array.
  f_full = (double*)calloc((L+1)*(2*L), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)
  for (t=1; t<=L-1; t++)
    memcpy(&f_full[t*f_stride], &f[(t-1)*f_stride], 
	   (2*L)*sizeof(double));

  // Define poles for all phi.
  for (p=0; p<=2*L-1; p++) {
    f_full[0*f_stride + p] = f_np; 
    f_full[L*f_stride + p] = f_sp; 
  }

  // Perform forward transform.
  ssht_adjoint_mw_inverse_sov_sym_ss_real(flm, f_full, L, 
					  dl_method, verbosity);

  // Free memory.
  free(f_full);

}
