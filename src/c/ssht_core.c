// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


/*!
 * \file ssht_core.c
 * Core algorithms to perform spin spherical harmonic transform on the sphere.
 *
 * \author Jason McEwen
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <tgmath.h>
#include <fftw3.h>

#include "ssht/ssht_types.h"
#include "ssht/ssht_error.h"
#include "ssht/ssht_dl.h"
#include "ssht/ssht_sampling.h"
#include "ssht/ssht_core.h"

#define MAX(a,b) ((a) > (b) ? (a) : (b))


//============================================================================
// MW algorithms
//============================================================================


/*!
 * Compute inverse transform for MW method using separation of
 * variables, fast Fourier transforms and exploiting all symmetries
 * (for complex spin signal).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_inverse_sov_sym(ssht_complex_double *f, const ssht_complex_double *flm,
				  int L, int spin,
				  ssht_dl_method_t dl_method,
				  int verbosity) {
    ssht_core_mw_lb_inverse_sov_sym(f, flm,
                                    0, L, spin,
                                    dl_method,
                                    verbosity);
}

/*!
 * Compute inverse transform for MW method using separation of
 * variables, fast Fourier transforms and exploiting all symmetries
 * (for complex spin signal).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_lb_inverse_sov_sym(ssht_complex_double *f, const ssht_complex_double *flm,
          int L0, int L, int spin,
          ssht_dl_method_t dl_method,
          int verbosity) {

  int el, m, mm;
  //int t, p;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, m_offset;
  double ssign, elfactor;
  ssht_complex_double mmfactor;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *exps, *m_factors;
  int exps_offset;
  double elmmsign, elssign;
  int spinneg;
  ssht_complex_double *Fmm, *fext;
  int Fmm_offset, Fmm_stride, fext_stride;
  fftw_plan plan;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  exps = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(exps)
  m_factors = calloc(2*L-1, sizeof *m_factors);
  SSHT_ERROR_MEM_ALLOC_CHECK(m_factors)

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
	   "Computing inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mw_inverse_sov_sym...");
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
  m_offset = L-1;
  for (el=MAX(L0, abs(spin)); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
    elfactor = ssign * sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el*el + el;
    for (m=-el; m<=el; m++)
      m_factors[m + m_offset] = flm[el2pel + m] * exps[m + exps_offset];

    for (mm=0; mm<=el; mm++) {
      double mm_factor;
      int mm_offset = mm*dl_stride;
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;
      mm_factor = elfactor * elssign * dl[mm_offset - spinneg + dl_offset];
      for (m=-el; m<=-1; m++) {
    	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset] +=
          mm_factor
          * m_factors[m + m_offset]
    	  * elmmsign * dl[mm_offset - m + dl_offset];
      }
      for (m=0; m<=el; m++) {
    	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset] +=
    	  mm_factor
          * m_factors[m + m_offset]
    	  * dl[mm_offset + m + dl_offset];
      }

    }

  }

  // Free dl memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);

  // Use symmetry to compute Fmm for negative mm.
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=-(L-1); m<=L-1; m++)
      Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset] =
	signs[abs(m)] * ssign
	* Fmm[(-mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];

  // Apply phase modulation to account for sampling offset.
  for (mm=-(L-1); mm<=L-1; mm++) {
    mmfactor = cexp(I*mm*SSHT_PI/(2.0*L-1.0));
    for (m=-(L-1); m<=L-1; m++)
      Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset] *=
  	mmfactor;
  }

  // Allocate space for function values.
  fext = (ssht_complex_double*)calloc((2*L-1)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(fext)
  fext_stride = 2*L-1;

  // Apply spatial shift.
  for (mm=0; mm<=L-1; mm++)
    for (m=0; m<=L-1; m++)
      fext[mm*fext_stride + m] =
	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];
  for (mm=0; mm<=L-1; mm++)
    for (m=-(L-1); m<=-1; m++)
      fext[mm*fext_stride + (m+2*L-1)] =
	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=0; m<=L-1; m++)
      fext[(mm + 2*L-1)*fext_stride + m] =
	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=-(L-1); m<=-1; m++)
      fext[(mm+2*L-1)*fext_stride + m + 2*L-1] =
	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];

  // Perform 2D FFT.
  plan = fftw_plan_dft_2d(2*L-1, 2*L-1, Fmm, Fmm,
			  FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute_dft(plan, fext, fext);
  fftw_destroy_plan(plan);

  // Free Fmm memory.
  free(Fmm);

  // Extract f from version of f extended to the torus (fext).
  memcpy(f, fext, L*(2*L-1)*sizeof(ssht_complex_double));
  /* Memcpy equivalent to:
  for (t=0; t<=L-1; t++)
    for (p=0; p<=2*L-2; p++)
      f[t*fext_stride + p] = fext[t*fext_stride + p];
  */

  // Free fext memory.
  free(fext);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(m_factors);

}

/*!
 * Compute inverse transform for MW method of real scalar signal using
 * separation of variables, fast Fourier transforms and exploiting all
 * symmetries (including additional symmetries for real signals).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_inverse_sov_sym_real(double *f, const ssht_complex_double *flm,
				       int L,
				       ssht_dl_method_t dl_method,
				       int verbosity) {
    ssht_core_mw_lb_inverse_sov_sym_real(f, flm,
                                         0, L,
                                         dl_method,
                                         verbosity);
}

/*!
 * Compute inverse transform for MW method of real scalar signal using
 * separation of variables, fast Fourier transforms and exploiting all
 * symmetries (including additional symmetries for real signals).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_lb_inverse_sov_sym_real(double *f, const ssht_complex_double *flm,
				       int L0, int L,
				       ssht_dl_method_t dl_method,
				       int verbosity) {

  int el, m, mm;
  //int t, p;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, m_offset;
  double ssign, elfactor;
  ssht_complex_double *m_factors;
  ssht_complex_double mmfactor;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  ssht_complex_double *exps;
  int exps_offset;
  double elmmsign, elssign;
  int spinneg;
  ssht_complex_double *Fmm, *Fmm_shift;
  double *fext_real;
  int Fmm_offset, Fmm_stride;
	// int fext_stride;
  fftw_plan plan;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  exps = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(exps)
  m_factors = calloc(2*L-1, sizeof *m_factors);
  SSHT_ERROR_MEM_ALLOC_CHECK(m_factors)

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
	   "Computing inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mw_inverse_sov_sym_real...");
  }

  // Compute Fmm.
  Fmm = (ssht_complex_double*)calloc((2*L-1)*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_offset = L-1;
  Fmm_stride = L;
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
  m_offset = L-1;
  for (el=MAX(L0, abs(spin)); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
    elfactor = ssign * sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    el2pel = el*el + el;
    for (m=-el; m<=el; m++)
      m_factors[m + m_offset] = flm[el2pel + m] * exps[m + exps_offset];

    for (mm=0; mm<=el; mm++) {
      double mm_factor;
      int mm_offset = mm*dl_stride;
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      mm_factor = elfactor * elssign * dl[mm_offset - spinneg + dl_offset];

      for (m=0; m<=el; m++) {
        Fmm[(mm + Fmm_offset)*Fmm_stride + m] +=
          mm_factor
          * m_factors[m + m_offset]
          * dl[mm_offset + m + dl_offset];
      }

    }

  }

  // Free dl memory.
  free(dl);
  if (dl_method == SSHT_DL_RISBO)
    free(dl8);

  // Use symmetry to compute Fmm for negative mm.
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=0; m<=L-1; m++)
      Fmm[(mm + Fmm_offset)*Fmm_stride + m] =
	signs[abs(m)] * ssign
	* Fmm[(-mm + Fmm_offset)*Fmm_stride + m];

  // Apply phase modulation to account for sampling offset.
  for (mm=-(L-1); mm<=L-1; mm++) {
    mmfactor = cexp(I*mm*SSHT_PI/(2.0*L-1.0));
    for (m=0; m<=L-1; m++)
      Fmm[(mm + Fmm_offset)*Fmm_stride + m] *=
	mmfactor;
  }

  // Apply spatial shift.
  Fmm_shift = (ssht_complex_double*)calloc((2*L-1)*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_shift)
  for (mm=0; mm<=L-1; mm++)
    for (m=0; m<=L-1; m++)
      Fmm_shift[mm*Fmm_stride + m] =
	Fmm[(mm + Fmm_offset)*Fmm_stride + m];
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=0; m<=L-1; m++)
      Fmm_shift[(mm + 2*L-1)*Fmm_stride + m] =
	Fmm[(mm + Fmm_offset)*Fmm_stride + m];

  // Allocate space for function values.
  fext_real = (double*)calloc((2*L-1)*(2*L-1), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(fext_real)
  // fext_stride = 2*L-1;

  // Perform 2D FFT.
  plan = fftw_plan_dft_c2r_2d(2*L-1, 2*L-1, Fmm_shift, fext_real,
			      FFTW_ESTIMATE);
  fftw_execute_dft_c2r(plan, Fmm_shift, fext_real);
  fftw_destroy_plan(plan);

  // Free Fmm memory.
  free(Fmm);
  free(Fmm_shift);

  // Extract f from version of f extended to the torus (fext).
  memcpy(f, fext_real, L*(2*L-1)*sizeof(double));
  /* Memcpy equivalent to:
  for (t=0; t<=L-1; t++)
    for (p=0; p<=2*L-2; p++)
      f[t*fext_stride + p] = fext[t*fext_stride + p];
  */

  // Free fext memory.
  free(fext_real);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(m_factors);

}


/*!
 * Compute inverse transform using direct method for MW sampling.
 *
 * \warning This algorithm is very slow and is included for
 * verification purposes only.
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mwdirect_inverse(ssht_complex_double *f, const ssht_complex_double *flm,
				 int L, int spin, int verbosity) {

  int t, p, m, el, ind, eltmp;
  double *dl;
  double *sqrt_tbl;
  double theta, phi, elfactor;
  int ssign;
  int dl_offset, dl_stride, f_stride;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)

    // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  ssign = spin & 1;
  ssign = 1 - ssign - ssign; // (-1)^spin

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mwdirect_inverse...");
  }

  // Initialise f with zeros.
  f_stride = 2*L-1;
  for (t=0; t<=L-1; t++)
    for (p=0; p<=2*L-2; p++)
      f[t*f_stride + p] = 0.0;

  // Compute inverse transform.
  dl = ssht_dl_calloc(L, SSHT_DL_FULL);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_FULL);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_FULL);
  for (t=0; t<=L-1; t++) {
    theta = ssht_sampling_mw_t2theta(t, L);
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
      if (el!=0 && el==abs(spin)) {
	for(eltmp=0; eltmp<=abs(spin); eltmp++)
	  ssht_dl_beta_risbo_full_table(dl, theta, L,
					SSHT_DL_FULL,
					eltmp, sqrt_tbl);
      }
      else {
	ssht_dl_beta_risbo_full_table(dl, theta, L,
				      SSHT_DL_FULL,
				      el, sqrt_tbl);
      }

      for (m=-el; m<=el; m++) {
	ssht_sampling_elm2ind(&ind, el, m);
	for (p=0; p<=2*L-2; p++) {
	  phi = ssht_sampling_mw_p2phi(p, L);
	  f[t*f_stride + p] +=
	    ssign
	    * elfactor
	    * cexp(I*m*phi)
	    * dl[(m+dl_offset)*dl_stride - spin + dl_offset]
	    * flm[ind];
	}
      }

    }
  }

  free(sqrt_tbl);
  free(dl);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

}


/*!
 * Compute inverse transform using direct method with separation of
 * variables for MW sampling.
 *
 * \warning This algorithm is included for verification purposes only.
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mwdirect_inverse_sov(ssht_complex_double *f, const ssht_complex_double *flm,
				     int L, int spin, int verbosity) {

  int t, p, m, el, ind;
  int ftm_stride, ftm_offset, f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  double *sqrt_tbl, *signs;
  ssht_complex_double *ftm, *inout;
  double theta, ssign, elfactor;
  fftw_plan plan;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)

  // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  for (m=0; m<=L-1; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }
  ssign = signs[abs(spin)];

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mwdirect_inverse_sov...");
  }

  // Compute ftm.
  ftm = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(ftm)
  ftm_stride = 2*L-1;
  ftm_offset = L-1;
  dlm1p1_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
  dl_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  for (t=0; t<=L-1; t++) {
    theta = ssht_sampling_mw_t2theta(t, L);
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));

      // Compute dl line for given spin.
      ssht_dl_beta_kostelec_line_table(dlm1p1_line, dl_line,
      				       theta, L, -spin, el,
      				       sqrt_tbl, signs);
      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

      for (m=-el; m<=el; m++) {
	ssht_sampling_elm2ind(&ind, el, m);
	ftm[t*ftm_stride + m + ftm_offset] +=
	  ssign
	  * elfactor
	  * dl_line[m + L-1]
	  * flm[ind];
      }
    }
  }

  // Free dl memory.
  free(dlm1p1_line);
  free(dl_line);

  // Compute f.
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  f_stride = 2*L-1;
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  for (t=0; t<=L-1; t++) {
    for (m=0; m<=L-1; m++)
      inout[m] = ftm[t*ftm_stride + m + ftm_offset];
    for (m=-(L-1); m<=-1; m++)
      inout[m+2*L-1] = ftm[t*ftm_stride + m + ftm_offset];
    fftw_execute_dft(plan, inout, inout);
    for (p=0; p<=2*L-2; p++)
      f[t*f_stride + p] = inout[p];
  }
  fftw_destroy_plan(plan);

  // Free memory.
  free(ftm);
  free(inout);
  free(signs);
  free(sqrt_tbl);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

}


/*!
 * Compute forward transform for MW method using separation of
 * variables, fast Fourier transforms, performing convolution with
 * weights as product in transformed space and exploiting all
 * symmetries (for complex spin signal).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_forward_sov_conv_sym(ssht_complex_double *flm, const ssht_complex_double *f,
                       int L, int spin,
                       ssht_dl_method_t dl_method,
                       int verbosity) {
    ssht_core_mw_lb_forward_sov_conv_sym(flm, f,
                                         0, L, spin,
                                         dl_method,
                                         verbosity);
}

/*!
 * Compute forward transform for MW method using separation of
 * variables, fast Fourier transforms, performing convolution with
 * weights as product in transformed space and exploiting all
 * symmetries (for complex spin signal).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_lb_forward_sov_conv_sym(ssht_complex_double *flm, const ssht_complex_double *f,
				       int L0, int L, int spin,
				       ssht_dl_method_t dl_method,
				       int verbosity) {

  int el, m, mm, ind, t, r;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel;
  double ssign, elfactor;
  fftw_plan plan, plan_bwd, plan_fwd;
  ssht_complex_double *inout;
  ssht_complex_double *Fmt, *Fmm, *Gmm, *m_mm_factor;
  ssht_complex_double *w, *wr;
  ssht_complex_double *Fmm_pad, *tmp_pad;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  int w_offset;
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
	   "Computing forward transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mw_forward_sov_conv_sym...");
  }

  // Compute Fourier transform over phi, i.e. compute Fmt.
  Fmt = (ssht_complex_double*)calloc((2*L-1)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmt)
  Fmt_stride = 2*L-1;
  Fmt_offset = L-1;
  f_stride = 2*L-1;
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_ESTIMATE);
  for (t=0; t<=L-1; t++) {
    memcpy(inout, &f[t*f_stride], f_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(m=0; m<=L-1; m++)
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m] / (2.0*L-1.0);
    for(m=-(L-1); m<=-1; m++)
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m+2*L-1] / (2.0*L-1.0);
  }

  // Extend Fmt periodically.
  for (m=-(L-1); m<=L-1; m++)
    for (t=L; t<=2*L-2; t++)
      Fmt[(m+Fmt_offset)*Fmt_stride + t] =
	signs[abs(m)] * ssign * Fmt[(m+Fmt_offset)*Fmt_stride + (2*L-2-t)];

  // Compute Fourier transform over theta, i.e. compute Fmm.
  Fmm = (ssht_complex_double*)calloc((2*L-1)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L-1;
  Fmm_offset = L-1;
  for (m=-(L-1); m<=L-1; m++) {
    memcpy(inout, &Fmt[(m+Fmt_offset)*Fmt_stride], Fmt_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L-1; mm++)
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] =
	inout[mm] / (2.0*L-1.0);
    for(mm=-(L-1); mm<=-1; mm++)
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] =
	inout[mm+2*L-1] / (2.0*L-1.0);
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (m=-(L-1); m<=L-1; m++)
    for (mm=-(L-1); mm<=L-1; mm++)
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] *=
  	expsmm[mm + exps_offset];

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
      Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] =
	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);

  }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);

  // Precompute factors depending on particular Gmm, to be used later
  // in computing the flm.
  m_mm_factor = calloc(L*(2*L-1), sizeof *m_mm_factor);
  SSHT_ERROR_MEM_ALLOC_CHECK(m_mm_factor)

  for (mm = 1; mm < L; mm++) {
    for (m = -L+1; m < L; m++) {
      m_mm_factor[mm*Fmm_stride + m + Fmm_offset] =
        expsm[m + exps_offset]
        * ( Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]
            + signs[abs(m)] * ssign
              * Gmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] );
    }
  }

  // Compute flm.
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
  for (el=0; el<=L-1; el++) {
    for (m=-el; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (el=MAX(L0, abs(spin)); el<=L-1; el++) {
    double mm_factor;

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
    elssign = spin <= 0 ? 1.0 : signs[el];

    mm_factor = ssign * elfactor * elssign * dl[- spinneg + dl_offset];

    for (m=-el; m<=-1; m++) {
      // mm = 0
      flm[m + el2pel] +=
        mm_factor
        * signs[el]
        * expsm[m + exps_offset]
        * dl[0*dl_stride - m + dl_offset]
        * Gmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
    }
    for (m=0; m<=el; m++) {
      // mm = 0
      flm[m + el2pel] +=
        mm_factor
        * expsm[m + exps_offset]
        * dl[0*dl_stride + m + dl_offset]
        * Gmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
    }

    for (mm=1; mm<=el; mm++) {
      int mm_offset = mm * dl_stride;
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      mm_factor =
          ssign
          * elfactor
          * elssign
          * dl[mm_offset - spinneg + dl_offset];

      for (m=-el; m<=-1; m++) {
        flm[m + el2pel] +=
          mm_factor
          * elmmsign * dl[mm_offset - m + dl_offset]
          * m_mm_factor[mm*Fmm_stride + m + Fmm_offset];
      }
      for (m=0; m<=el; m++) {
        flm[m + el2pel] +=
          mm_factor
          * dl[mm_offset + m + dl_offset]
          * m_mm_factor[mm*Fmm_stride + m + Fmm_offset];
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
  free(m_mm_factor);
  free(sqrt_tbl);
  free(signs);
  free(expsm);
  free(expsmm);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Forward transform computed!");

}


/*!
 * Compute forward transform for MW method of real scalar signal using
 * separation of variables, fast Fourier transforms, performing
 * convolution with weights as product in transformed space and
 * exploiting all symmetries (including additional symmetries for real
 * signals).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_forward_sov_conv_sym_real(ssht_complex_double *flm, const double *f,
                        int L,
                        ssht_dl_method_t dl_method,
                        int verbosity) {
    ssht_core_mw_lb_forward_sov_conv_sym_real(flm, f,
                                              0, L,
                                              dl_method,
                                              verbosity);
}

/*!
 * Compute forward transform for MW method of real scalar signal using
 * separation of variables, fast Fourier transforms, performing
 * convolution with weights as product in transformed space and
 * exploiting all symmetries (including additional symmetries for real
 * signals).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_lb_forward_sov_conv_sym_real(ssht_complex_double *flm, const double *f,
					    int L0, int L,
					    ssht_dl_method_t dl_method,
					    int verbosity) {

  int el, m, mm, ind, ind_nm, t, r;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel;
  double ssign, elfactor;
  fftw_plan plan, plan_bwd, plan_fwd;
  double *in_real;
  ssht_complex_double *inout, *out;
  ssht_complex_double *Fmt, *Fmm, *Gmm, *m_mm_factor;
  ssht_complex_double *w, *wr;
  ssht_complex_double *Fmm_pad, *tmp_pad;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset, Gmm_stride;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  int w_offset;
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
	   "Computing forward transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mw_forward_sov_conv_sym_real...");
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
      Fmt[m*Fmt_stride + t] = out[m] / (2.0*L-1.0);
  }
  free(in_real);
  free(out);
  fftw_destroy_plan(plan);

  // Extend Fmt periodically.
  for (m=0; m<=L-1; m++)
    for (t=L; t<=2*L-2; t++)
      Fmt[m*Fmt_stride + t] =
	signs[abs(m)] * ssign * Fmt[m*Fmt_stride + (2*L-2-t)];

  // Compute Fourier transform over theta, i.e. compute Fmm.
  Fmm = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L-1;
  Fmm_offset = L-1;
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (m=0; m<=L-1; m++) {
    memcpy(inout, &Fmt[m*Fmt_stride], Fmt_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L-1; mm++)
      Fmm[m*Fmm_stride + mm + Fmm_offset] =
	inout[mm] / (2.0*L-1.0);
    for(mm=-(L-1); mm<=-1; mm++)
      Fmm[m*Fmm_stride + mm + Fmm_offset] =
	inout[mm+2*L-1] / (2.0*L-1.0);
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Apply phase modulation to account for sampling offset.
  for (m=0; m<=L-1; m++)
    for (mm=-(L-1); mm<=L-1; mm++)
      Fmm[m*Fmm_stride + mm + Fmm_offset] *=
	expsmm[mm + exps_offset];

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
  Gmm = (ssht_complex_double*)calloc((2*L-1)*L, sizeof(ssht_complex_double));
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

  // Precompute factors depending on particular Gmm, to be used later
  // in computing the flm.
  m_mm_factor = calloc(L*L, sizeof *m_mm_factor);
  SSHT_ERROR_MEM_ALLOC_CHECK(m_mm_factor)

  for (mm = 1; mm < L; mm++) {
    for (m = 0; m < L; m++) {
      m_mm_factor[mm*Gmm_stride + m] =
        expsm[m]
        * ( Gmm[(mm+Fmm_offset)*Gmm_stride + m]
            + signs[m] * ssign
              * Gmm[(-mm+Fmm_offset)*Gmm_stride + m] );
    }
  }

  // Compute flm.
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
  for (el=0; el<=L-1; el++) {
    for (m=0; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (el=MAX(L0, abs(spin)); el<=L-1; el++) {
    double mm_factor;

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
    elssign = spin <= 0 ? 1.0 : signs[el];

    mm_factor = ssign * elfactor * elssign * dl[- spinneg + dl_offset];

    for (m=0; m<=el; m++) {
      // mm = 0
      flm[el2pel + m] +=
        mm_factor
        * expsm[m]
        * dl[0*dl_stride + m + dl_offset]
        * Gmm[(0+Fmm_offset)*Gmm_stride + m];
    }

    for (mm=1; mm<=el; mm++) {
      int mm_offset = mm * dl_stride;
      elmmsign = signs[el] * signs[mm];
      elssign = spin <= 0 ? 1.0 : elmmsign;

      mm_factor =
        ssign
        * elfactor
        * elssign
        * dl[mm_offset - spinneg + dl_offset];

      for (m=0; m<=el; m++) {
        flm[el2pel + m] +=
          mm_factor
          * dl[mm_offset + m + dl_offset]
          * m_mm_factor[mm*Gmm_stride + m];
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
  free(m_mm_factor);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Forward transform computed!");

}


//============================================================================
// MW South pole interfaces
//============================================================================


/*!
 * South pole wrapper for inverse transform for MW method.  The South
 * pole is defined by a single sample and its corresponding phi angle,
 * rather than specifying samples for all phi at the South pole (which
 * are simply related by the rotation of a spin function in its
 * tangent plane).
 *
 * \param[out] f Function on sphere (excluding South pole).
 * \param[out] f_sp Function sample on South pole.
 * \param[out] phi_sp Phi angle corresponding to quoted sample at
 * South pole.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_inverse_sov_sym_pole(ssht_complex_double *f,
				       ssht_complex_double *f_sp, double *phi_sp,
				       const ssht_complex_double *flm,
				       int L, int spin,
				       ssht_dl_method_t dl_method,
				       int verbosity) {

  ssht_complex_double* f_full;
  int f_stride = 2*L-1;

  // Allocate full array.
  f_full = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)

  // Perform inverse transform.
  ssht_core_mw_inverse_sov_sym(f_full, flm, L, spin,
			       dl_method, verbosity);

  // Copy output function values, including separate point for South pole.
  memcpy(f, f_full, (L-1)*(2*L-1)*sizeof(ssht_complex_double));
  *f_sp = f_full[(L-1)*f_stride + 0];
  *phi_sp = ssht_sampling_mw_p2phi(0, L);

  // Free memory.
  free(f_full);

}


/*!
 * South pole wrapper for inverse transform of real scalar function
 * for MW method.  The South pole is defined by a single sample,
 * rather than specifying samples for all phi at the South
 * pole (which for a scalar function are identical).
 *
 * \param[out] f Function on sphere (excluding South pole).
 * \param[out] f_sp Function sample on South pole.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_inverse_sov_sym_real_pole(double *f,
					    double *f_sp,
					    const ssht_complex_double *flm,
					    int L,
					    ssht_dl_method_t dl_method,
					    int verbosity) {

  double *f_full;
  int f_stride = 2*L-1;

  // Allocate full array.
  f_full = (double*)calloc(L*(2*L-1), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)

  // Perform inverse transform.
  ssht_core_mw_inverse_sov_sym_real(f_full, flm, L,
				    dl_method, verbosity);

  // Copy output function values, including separate point for South pole.
  memcpy(f, f_full, (L-1)*(2*L-1)*sizeof(double));
  *f_sp = f_full[(L-1)*f_stride + 0];

  // Free memory.
  free(f_full);

}


/*!
 * South pole wrapper for forward transform for MW method.  The South
 * pole is defined by a single sample and its corresponding phi angle,
 * rather than specifying samples for all phi at the South pole (which
 * are simply related by the rotation of a spin function in its
 * tangent plane).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere (excluding South pole).
 * \param[in] f_sp Function sample on South pole.
 * \param[in] phi_sp Phi angle corresponding to quoted sample at
 * South pole.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_forward_sov_conv_sym_pole(ssht_complex_double *flm, const ssht_complex_double *f,
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
  ssht_core_mw_forward_sov_conv_sym(flm, f_full, L, spin,
				    dl_method, verbosity);

  // Free memory.
  free(f_full);

}


/*!
 * South pole wrapper for forward transform of real scalar function
 * for MW method.  The South pole is defined by a single sample,
 * rather than specifying samples for all phi at the South
 * pole (which for a scalar function are identical).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere (excluding South pole).
 * \param[in] f_sp Function sample on South pole.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_forward_sov_conv_sym_real_pole(ssht_complex_double *flm,
						 const double *f,
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
  ssht_core_mw_forward_sov_conv_sym_real(flm, f_full, L,
					 dl_method, verbosity);

  // Free memory.
  free(f_full);

}


//============================================================================
// MW SS algorithms
//============================================================================

/*!
 * Compute inverse transform for MW method with symmetric sampling
 * using separation of variables, fast Fourier transforms and
 * exploiting all symmetries (for complex spin signal).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_inverse_sov_sym_ss(ssht_complex_double *f, const ssht_complex_double *flm,
                     int L, int spin,
                     ssht_dl_method_t dl_method,
                     int verbosity) {
    ssht_core_mw_lb_inverse_sov_sym_ss(f, flm,
                                       0, L, spin,
                                       dl_method,
                                       verbosity);
}

/*!
 * Compute inverse transform for MW method with symmetric sampling
 * using separation of variables, fast Fourier transforms and
 * exploiting all symmetries (for complex spin signal).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_lb_inverse_sov_sym_ss(ssht_complex_double *f, const ssht_complex_double *flm,
				     int L0, int L, int spin,
				     ssht_dl_method_t dl_method,
				     int verbosity) {

  int el, m, mm, ind;
  //int t, p;
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
  ssht_complex_double *Fmm, *fext;
  int Fmm_offset, Fmm_stride, fext_stride;
  fftw_plan plan;

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
	   "Computing inverse transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mw_inverse_sov_sym_ss...");
  }

  // Compute Fmm.
  // Note that m and mm indices are increased in size by one and
  // will be filled with zeros by calloc.
  Fmm = (ssht_complex_double*)calloc((2*L)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_offset = L-1;
  Fmm_stride = 2*L;
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
  inds_offset = L-1;
  for (el=MAX(L0, abs(spin)); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
    	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset] +=
    	  ssign
    	  * elfactor
	  * exps[m + exps_offset]
    	  * elmmsign * dl[mm*dl_stride - m + dl_offset]
    	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
    	  * flm[ind];
      }
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
    	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset] +=
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
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=-(L-1); m<=L; m++)
      Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset] =
	signs[abs(m)] * ssign
	* Fmm[(-mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];

  // Allocate space for function values.
  // Note that t and p indices of fext are increased in size by
  // one compared to usual sampling.
  fext = (ssht_complex_double*)calloc((2*L)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(fext)
  fext_stride = 2*L;

  // Apply spatial shift.
  for (mm=0; mm<=L; mm++)
    for (m=0; m<=L; m++)
      fext[mm*fext_stride + m] =
  	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];
  for (mm=0; mm<=L; mm++)
    for (m=-(L-1); m<=-1; m++)
      fext[mm*fext_stride + (m+2*L-1+1)] =
  	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=0; m<=L; m++)
      fext[(mm + 2*L-1+1)*fext_stride + m] =
  	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=-(L-1); m<=-1; m++)
      fext[(mm+2*L-1+1)*fext_stride + m + 2*L-1+1] =
  	Fmm[(mm + Fmm_offset)*Fmm_stride + m + Fmm_offset];

  // Perform 2D FFT.
  plan = fftw_plan_dft_2d(2*L, 2*L, Fmm, Fmm,
			  FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute_dft(plan, fext, fext);
  fftw_destroy_plan(plan);

  // Free Fmm memory.
  free(Fmm);

  // Extract f from version of f extended to the torus (fext).
  // Note that t and p indices of fext are increased in size by
  // one compared to usual sampling.
  memcpy(f, fext, (L+1)*(2*L)*sizeof(ssht_complex_double));

  // Free fext memory.
  free(fext);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(inds);

}

/*!
 * Compute inverse transform for MW method with symmetric sampling of
 * real scalar signal using separation of variables, fast Fourier
 * transforms and exploiting all symmetries (including additional
 * symmetries for real signals).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_inverse_sov_sym_ss_real(double *f, const ssht_complex_double *flm,
                      int L,
                      ssht_dl_method_t dl_method,
                      int verbosity) {
    ssht_core_mw_lb_inverse_sov_sym_ss_real(f, flm,
                                            0, L,
                                            dl_method,
                                            verbosity);
}

/*!
 * Compute inverse transform for MW method with symmetric sampling of
 * real scalar signal using separation of variables, fast Fourier
 * transforms and exploiting all symmetries (including additional
 * symmetries for real signals).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_lb_inverse_sov_sym_ss_real(double *f, const ssht_complex_double *flm,
					  int L0, int L,
					  ssht_dl_method_t dl_method,
					  int verbosity) {

  int el, m, mm, ind;
  //int t, p;
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
  ssht_complex_double *Fmm, *Fmm_shift;
  double *fext_real;
  int Fmm_offset, Fmm_stride;
  fftw_plan plan;
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
	   "Computing inverse transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mw_inverse_sov_sym_ss_real...");
  }

  // Compute Fmm.
  // Note that m and mm indices are increased in size by one and
  // will be filled with zeros by calloc.
  Fmm = (ssht_complex_double*)calloc((2*L)*(L+1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_offset = L-1;
  Fmm_stride = L+1;
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  if (dl_method == SSHT_DL_RISBO) {
    dl8 = ssht_dl_calloc(L, SSHT_DL_QUARTER_EXTENDED);
    SSHT_ERROR_MEM_ALLOC_CHECK(dl8)
  }
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);
  inds_offset = L-1;
  for (el=MAX(L0, abs(spin)); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
    	Fmm[(mm + Fmm_offset)*Fmm_stride + m] +=
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
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=0; m<=L; m++)
      Fmm[(mm + Fmm_offset)*Fmm_stride + m] =
	signs[abs(m)] * ssign
	* Fmm[(-mm + Fmm_offset)*Fmm_stride + m];

  // Apply spatial shift.
  Fmm_shift = (ssht_complex_double*)calloc((2*L)*(L+1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm_shift)
  for (mm=0; mm<=L-1; mm++)
    for (m=0; m<=L-1; m++)
      Fmm_shift[mm*Fmm_stride + m] =
	Fmm[(mm + Fmm_offset)*Fmm_stride + m];
  for (mm=-(L-1); mm<=-1; mm++)
    for (m=0; m<=L-1; m++)
      Fmm_shift[(mm + 2*L-1+1)*Fmm_stride + m] =
	Fmm[(mm + Fmm_offset)*Fmm_stride + m];

  // Allocate space for function values.
  // Note that t and p indices of fext are increased in size by
  // one compared to usual sampling.
  fext_real = (double*)calloc((2*L)*(2*L), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(fext_real)  

  // Perform 2D FFT.
  plan = fftw_plan_dft_c2r_2d(2*L, 2*L, Fmm_shift, fext_real,
			      FFTW_ESTIMATE);
  fftw_execute_dft_c2r(plan, Fmm_shift, fext_real);
  fftw_destroy_plan(plan);

  // Free Fmm memory.
  free(Fmm);
  free(Fmm_shift);

  // Extract f from version of f extended to the torus (fext).
  // Note that t and p indices of fext are increased in size by
  // one compared to usual sampling.
  memcpy(f, fext_real, (L+1)*(2*L)*sizeof(double));

  // Free fext memory.
  free(fext_real);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs);
  free(exps);
  free(inds);

}


/*!
 * Compute inverse transform using direct method for MW symmetric
 * sampling.
 *
 * \warning This algorithm is very slow and is included for
 * verification purposes only.
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mwdirect_inverse_ss(ssht_complex_double *f, const ssht_complex_double *flm,
				   int L, int spin, int verbosity) {

  int t, p, m, el, ind, eltmp;
  double *dl;
  double *sqrt_tbl;
  double theta, phi, elfactor;
  int ssign;
  int dl_offset, dl_stride, f_stride;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)

    // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  ssign = spin & 1;
  ssign = 1 - ssign - ssign; // (-1)^spin

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mwdirect_inverse_ss...");
  }

  // Initialise f with zeros.
  f_stride = 2*L;
  for (t=0; t<=L; t++)
    for (p=0; p<=2*L-1; p++)
      f[t*f_stride + p] = 0.0;

  // Compute inverse transform.
  dl = ssht_dl_calloc(L, SSHT_DL_FULL);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_FULL);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_FULL);
  for (t=0; t<=L; t++) {
    theta = ssht_sampling_mw_ss_t2theta(t, L);
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
      if (el!=0 && el==abs(spin)) {
	for(eltmp=0; eltmp<=abs(spin); eltmp++)
	  ssht_dl_beta_risbo_full_table(dl, theta, L,
					SSHT_DL_FULL,
					eltmp, sqrt_tbl);
      }
      else {
	ssht_dl_beta_risbo_full_table(dl, theta, L,
				      SSHT_DL_FULL,
				      el, sqrt_tbl);
      }

      for (m=-el; m<=el; m++) {
	ssht_sampling_elm2ind(&ind, el, m);
	for (p=0; p<=2*L-1; p++) {
	  phi = ssht_sampling_mw_ss_p2phi(p, L);
	  f[t*f_stride + p] +=
	    ssign
	    * elfactor
	    * cexp(I*m*phi)
	    * dl[(m+dl_offset)*dl_stride - spin + dl_offset]
	    * flm[ind];
	}
      }

    }
  }

  free(sqrt_tbl);
  free(dl);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

}

/*!
 * Compute forward transform for MW method with symmetric sampling
 * using separation of variables, fast Fourier transforms, performing
 * convolution with weights as product in transformed space and
 * exploiting all symmetries (for complex spin signal).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_forward_sov_conv_sym_ss(ssht_complex_double *flm, const ssht_complex_double *f,
                      int L, int spin,
                      ssht_dl_method_t dl_method,
                      int verbosity) {
    ssht_core_mw_lb_forward_sov_conv_sym_ss(flm, f,
                                            0, L, spin,
                                            dl_method,
                                            verbosity);
}

/*!
 * Compute forward transform for MW method with symmetric sampling
 * using separation of variables, fast Fourier transforms, performing
 * convolution with weights as product in transformed space and
 * exploiting all symmetries (for complex spin signal).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_lb_forward_sov_conv_sym_ss(ssht_complex_double *flm, const ssht_complex_double *f,
					  int L0, int L, int spin,
					  ssht_dl_method_t dl_method,
					  int verbosity) {

  int el, m, mm, ind, t, r;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan, plan_bwd, plan_fwd;
  ssht_complex_double *inout;
  ssht_complex_double *Fmt, *Fmm, *Gmm;
  ssht_complex_double *w, *wr;
  ssht_complex_double *Fmm_pad, *tmp_pad;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  int w_offset;
  ssht_complex_double *expsm, *expsmm;
  int exps_offset;
  int elmmsign, elssign;
  int spinneg;
  int Gmm_stride, Gmm_offset;

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
	   "Computing forward transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mw_forward_sov_conv_sym_ss...");
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
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m] / (2.0*L);
    for(m=-(L-1); m<=-1; m++)
      Fmt[(m+Fmt_offset)*Fmt_stride + t] = inout[m+2*L-1+1] / (2.0*L);
  }
  fftw_destroy_plan(plan);
  free(inout);

  // Extend Fmt periodically.
  for (m=-(L-1); m<=L; m++)
    for (t=L+1; t<=2*L-1; t++)
      Fmt[(m+Fmt_offset)*Fmt_stride + t] =
	signs[abs(m)] * ssign * Fmt[(m+Fmt_offset)*Fmt_stride + (2*L-t)];

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
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] =
	inout[mm] / (2.0*L);
    for(mm=-(L-1); mm<=-1; mm++)
      Fmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] =
	inout[mm+2*L-1+1] / (2.0*L);
  }
  fftw_destroy_plan(plan);
  free(inout);

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
  for (el=MAX(L0, abs(spin)); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	* Gmm[(0+Gmm_offset)*Gmm_stride + m + Gmm_offset];
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
	* Gmm[(0+Gmm_offset)*Gmm_stride + m + Gmm_offset];
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
      	  * ( Gmm[(mm+Gmm_offset)*Gmm_stride + m + Gmm_offset]
      	      + signs[-m] * ssign
      	      * Gmm[(-mm+Gmm_offset)*Gmm_stride + m + Gmm_offset]);
      }
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign
	  * elfactor
	  * expsm[m + exps_offset]
	  * dl[mm*dl_stride + m + dl_offset]
	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
	  * ( Gmm[(mm+Gmm_offset)*Gmm_stride + m + Gmm_offset]
	      + signs[m] * ssign
	      * Gmm[(-mm+Gmm_offset)*Gmm_stride + m + Gmm_offset]);
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
    printf("%s%s", SSHT_PROMPT, "Forward transform computed!");

}

/*!
 * Compute forward transform for MW method using symmetric sampling of
 * real scalar signal using separation of variables, fast Fourier
 * transforms, performing convolution with weights as product in
 * transformed space and exploiting all symmetries (including
 * additional symmetries for real signals).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_forward_sov_conv_sym_ss_real(ssht_complex_double *flm, const double *f,
                           int L,
                           ssht_dl_method_t dl_method,
                           int verbosity) {
    ssht_core_mw_lb_forward_sov_conv_sym_ss_real(flm, f,
                                                 0, L,
                                                 dl_method,
                                                 verbosity);
}

/*!
 * Compute forward transform for MW method using symmetric sampling of
 * real scalar signal using separation of variables, fast Fourier
 * transforms, performing convolution with weights as product in
 * transformed space and exploiting all symmetries (including
 * additional symmetries for real signals).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_lb_forward_sov_conv_sym_ss_real(ssht_complex_double *flm, const double *f,
					       int L0, int L,
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
  ssht_complex_double *inout, *out;
  ssht_complex_double *Fmt, *Fmm, *Gmm;
  ssht_complex_double *w, *wr;
  ssht_complex_double *Fmm_pad, *tmp_pad;
  int f_stride, Fmt_stride, Fmt_offset, Fmm_stride, Fmm_offset, Gmm_stride;
  double *dl;
  double *dl8 = NULL;
  int dl_offset, dl_stride;
  int w_offset;
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
	   "Computing forward transform using MW symmetric sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_mw_forward_sov_conv_sym_ss_real...");
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
      Fmt[m*Fmt_stride + t] = out[m] / (2.0*L);

  }
  free(in_real);
  free(out);
  fftw_destroy_plan(plan);

  // Extend Fmt periodically.
  for (m=0; m<=L; m++)
    for (t=L+1; t<=2*L-1; t++)
      Fmt[m*Fmt_stride + t] =
	signs[abs(m)] * ssign * Fmt[m*Fmt_stride + (2*L-t)];

  // Compute Fourier transform over theta, i.e. compute Fmm.
  // Note that m and mm indices are increased in size by one.
  Fmm = (ssht_complex_double*)calloc((L+1)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L;
  Fmm_offset = L-1;
  inout = (ssht_complex_double*)calloc(2*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (m=0; m<=L; m++) {
    memcpy(inout, &Fmt[m*Fmt_stride], Fmt_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(mm=0; mm<=L; mm++)
      Fmm[m*Fmm_stride + mm + Fmm_offset] =
	inout[mm] / (2.0*L);
    for(mm=-(L-1); mm<=-1; mm++)
      Fmm[m*Fmm_stride + mm + Fmm_offset] =
	inout[mm+2*L-1+1] / (2.0*L);
  }
  fftw_destroy_plan(plan);
  free(inout);

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
  Gmm = (ssht_complex_double*)calloc((2*L-1)*L, sizeof(ssht_complex_double));
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
  for (el=MAX(L0, abs(spin)); el<=L-1; el++) {

    // Compute Wigner plane.
    switch (dl_method) {

      case SSHT_DL_RISBO:
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	if (el!=0 && el==MAX(L0, abs(spin))) {
	  for(eltmp=0; eltmp<=MAX(L0, abs(spin)); eltmp++)
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
	* Gmm[(0+Fmm_offset)*Gmm_stride + m];
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
	  * ( Gmm[(mm+Fmm_offset)*Gmm_stride + m]
	      + signs[m] * ssign
	      * Gmm[(-mm+Fmm_offset)*Gmm_stride + m]);
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
    printf("%s%s", SSHT_PROMPT, "Forward transform computed!");

}


//============================================================================
// MW SS Noth-South pole interfaces
//============================================================================


/*!
 * North-South pole wrapper for inverse transform for MW method with
 * symmetric sampling.  The poles are defined by single samples and
 * their corresponding phi angle, rather than specifying samples for
 * all phi at the poles (which are simply related by the rotation of a
 * spin function in its tangent plane).
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
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_inverse_sov_sym_ss_pole(ssht_complex_double *f,
					  ssht_complex_double *f_np, double *phi_np,
					  ssht_complex_double *f_sp, double *phi_sp,
					  const ssht_complex_double *flm,
					  int L, int spin,
					  ssht_dl_method_t dl_method,
					  int verbosity) {

  ssht_complex_double* f_full;
  int t, f_stride = 2*L;

  // Allocate full array.
  f_full = (ssht_complex_double*)calloc((L+1)*(2*L), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)

  // Perform inverse transform.
  ssht_core_mw_inverse_sov_sym_ss(f_full, flm, L, spin,
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
 * North-South pole wrapper for inverse transform of real scalar
 * function for MW method with symmetric sampling.  The poles are
 * defined by single samples, rather than specifying samples for all
 * phi at the poles (which for a scalar function are identical).
 *
 * \param[out] f Function on sphere (excluding poles).
 * \param[out] f_sp Function sample on South pole.
 * \param[out] f_np Function sample on North pole.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_inverse_sov_sym_ss_real_pole(double *f,
					       double *f_np,
					       double *f_sp,
					       const ssht_complex_double *flm,
					       int L,
					       ssht_dl_method_t dl_method,
					       int verbosity) {

  double *f_full;
  int t, f_stride = 2*L;

  // Allocate full array.
  f_full = (double*)calloc((L+1)*(2*L), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_full)

  // Perform inverse transform.
  ssht_core_mw_inverse_sov_sym_ss_real(f_full, flm, L,
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
 * North-South pole wrapper for forward transform for MW method with
 * symmetric sampling.  The poles are defined by single samples and their
 * corresponding phi angle, rather than specifying samples for all phi
 * at the poles (which are simply related by the rotation of a
 * spin function in its tangent plane).
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
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_forward_sov_conv_sym_ss_pole(ssht_complex_double *flm, const ssht_complex_double *f,
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
  ssht_core_mw_forward_sov_conv_sym_ss(flm, f_full, L, spin,
				       dl_method, verbosity);

  // Free memory.
  free(f_full);

}


/*!
 * North-South pole wrapper for forward transform of real scalar
 * function for MW method with symmetric sampling.  The poles are
 * defined by single samples, rather than specifying samples for all
 * phi at the poles (which for a scalar function are identical).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere (excluding poles).
 * \param[in] f_np Function sample on North pole.
 * \param[in] f_sp Function sample on South pole.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] dl_method Method to use when compute Wigner functions.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_mw_forward_sov_conv_sym_ss_real_pole(ssht_complex_double *flm,
						    const double *f,
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
  ssht_core_mw_forward_sov_conv_sym_ss_real(flm, f_full, L,
					    dl_method, verbosity);

  // Free memory.
  free(f_full);

}


//============================================================================
// GL algorithms
//============================================================================


/*!
 * Compute inverse transform using direct method with separation of
 * variables for GL sampling.
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_gl_inverse_sov(ssht_complex_double *f, const ssht_complex_double *flm,
			      int L, int spin, int verbosity) {

  int t, p, m, el, ind;
  int ftm_stride, ftm_offset, f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  double *sqrt_tbl, *signs;
  ssht_complex_double *ftm, *inout;
  double theta, ssign, elfactor;
  fftw_plan plan;
  double *thetas, *weights;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)

  // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  for (m=0; m<=L-1; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }
  ssign = signs[abs(spin)];

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing inverse transform using GL sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_gl_inverse_sov...");
  }

  // Compute weights and theta positions.
  thetas = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(thetas)
  weights = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(weights)
  ssht_sampling_gl_thetas_weights(thetas, weights, L);

  // Compute ftm.
  ftm = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(ftm)
  ftm_stride = 2*L-1;
  ftm_offset = L-1;
  dlm1p1_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
  dl_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  for (t=0; t<=L-1; t++) {
    theta = thetas[t];
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));

      // Compute dl line for given spin.
      ssht_dl_beta_kostelec_line_table(dlm1p1_line, dl_line,
      				       theta, L, -spin, el,
      				       sqrt_tbl, signs);
      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

      for (m=-el; m<=el; m++) {
	ssht_sampling_elm2ind(&ind, el, m);
	ftm[t*ftm_stride + m + ftm_offset] +=
	  ssign
	  * elfactor
	  * dl_line[m + L-1]
	  * flm[ind];
      }
    }
  }

  // Free dl memory.
  free(dlm1p1_line);
  free(dl_line);

  // Free memory for thetas and weights.
  free(thetas);
  free(weights);

  // Compute f.
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  f_stride = 2*L-1;
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  for (t=0; t<=L-1; t++) {
    for (m=0; m<=L-1; m++)
      inout[m] = ftm[t*ftm_stride + m + ftm_offset];
    for (m=-(L-1); m<=-1; m++)
      inout[m+2*L-1] = ftm[t*ftm_stride + m + ftm_offset];
    fftw_execute_dft(plan, inout, inout);
    for (p=0; p<=2*L-2; p++)
      f[t*f_stride + p] = inout[p];
  }
  fftw_destroy_plan(plan);

  // Free memory.
  free(ftm);
  free(inout);
  free(signs);
  free(sqrt_tbl);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

}


/*!
 * Compute inverse transform of real scalar signal using direct method
 * with separation of variables for GL sampling (symmetries for real
 * signals are exploited).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_gl_inverse_sov_real(double *f, const ssht_complex_double *flm,
				   int L, int verbosity) {

  int t, p, m, el, ind;
  int ftm_stride, ftm_offset, f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  double *sqrt_tbl, *signs;
  ssht_complex_double *ftm;
  ssht_complex_double *in;
  double *out;
  double theta, ssign, elfactor;
  fftw_plan plan;
  double *thetas, *weights;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)

  // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  for (m=0; m<=L-1; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }
  ssign = signs[abs(spin)];

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing inverse transform using GL sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_gl_inverse_sov_real...");
  }

  // Compute weights and theta positions.
  thetas = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(thetas)
  weights = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(weights)
  ssht_sampling_gl_thetas_weights(thetas, weights, L);

  // Compute ftm.
  ftm = (ssht_complex_double*)calloc(L*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(ftm)
  ftm_stride = L;
  ftm_offset = 0;
  dlm1p1_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
  dl_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  for (t=0; t<=L-1; t++) {
    theta = thetas[t];
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));

      // Compute half dl line for given spin.
      ssht_dl_beta_kostelec_halfline_table(dlm1p1_line, dl_line,
					   theta, L, -spin, el,
					   sqrt_tbl, signs);
      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

      for (m=0; m<=el; m++) {
	ssht_sampling_elm2ind(&ind, el, m);
	ftm[t*ftm_stride + m + ftm_offset] +=
	  ssign
	  * elfactor
	  * dl_line[m]
	  * flm[ind];
      }
    }
  }

  // Free dl memory.
  free(dlm1p1_line);
  free(dl_line);

  // Free memory for thetas and weights.
  free(thetas);
  free(weights);

  // Compute f.
  in = (ssht_complex_double*)calloc(L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in)
  out = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out)
  plan = fftw_plan_dft_c2r_1d(2*L-1, in, out, FFTW_MEASURE);
  f_stride = 2*L-1;
  for (t=0; t<=L-1; t++) {
    memcpy(in, &ftm[t*ftm_stride], L*sizeof(ssht_complex_double));
    fftw_execute_dft_c2r(plan, in, out);
    for (p=0; p<=2*L-2; p++)
      f[t*f_stride + p] = out[p];
  }
  fftw_destroy_plan(plan);

  // Free memory.
  free(ftm);
  free(in);
  free(out);
  free(signs);
  free(sqrt_tbl);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

}


/*!
 * Compute forward transform using Gauss-Legendgre quadrature with separation of
 * variables.
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_gl_forward_sov(ssht_complex_double *flm, const ssht_complex_double *f,
			      int L, int spin, int verbosity) {

  int t, m, el, ind;
  int f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  int el2pel, inds_offset;
  int *inds;
  double *sqrt_tbl, *signs;
  int Ftm_stride, Ftm_offset;
  ssht_complex_double *Ftm, *inout;
  double theta, ssign, elfactor;
  fftw_plan plan;
  double *thetas, *weights;
  double w;

   // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
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

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing forward transform using GL sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_gl_forward_sov...");
  }

  // Compute weights and theta positions.
  thetas = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(thetas)
  weights = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(weights)
  ssht_sampling_gl_thetas_weights(thetas, weights, L);

  // Compute Fourier transform over phi, i.e. compute Ftm.
  Ftm = (ssht_complex_double*)calloc(L*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = 2*L-1;
  Ftm_offset = L-1;
  f_stride = 2*L-1;
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (t=0; t<=L-1; t++) {
    memcpy(inout, &f[t*f_stride], f_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(m=0; m<=L-1; m++)
      Ftm[t*Ftm_stride + m + Ftm_offset] =
	inout[m] * 2.0 * SSHT_PI / (2.0*L-1.0) ;
    for(m=-(L-1); m<=-1; m++)
      Ftm[t*Ftm_stride + m + Ftm_offset] =
	inout[m+2*L-1] * 2.0 * SSHT_PI / (2.0*L-1.0);
  }
  fftw_destroy_plan(plan);

  // Compute flm.
  dlm1p1_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
  dl_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  inds_offset = L-1;
  for (el=0; el<=L-1; el++) {
    for (m=-el; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (t=0; t<=L-1; t++) {
    theta = thetas[t];
    w = weights[t];
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
      el2pel = el *el + el;
      for (m=-el; m<=el; m++)
	inds[m + inds_offset] = el2pel + m;

      // Compute dl line for given spin.
      ssht_dl_beta_kostelec_line_table(dlm1p1_line, dl_line,
      				       theta, L, -spin, el,
      				       sqrt_tbl, signs);
      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

      for (m=-el; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign
	  * elfactor
	  * w
	  * dl_line[m+L-1]
	  * Ftm[t*Ftm_stride + m + Ftm_offset];
      }
    }
  }

  // Free memory.
  free(dlm1p1_line);
  free(dl_line);
  free(thetas);
  free(weights);
  free(Ftm);
  free(inout);
  free(signs);
  free(sqrt_tbl);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Forward transform computed!");

}


/*!
 * Compute forward transform of real scalar signal using
 * Gauss-Legendgre quadrature with separation of variables (symmetries
 * for real signals are exploited).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_gl_forward_sov_real(ssht_complex_double *flm, const double *f,
				   int L, int verbosity) {

  int t, m, el, ind, ind_nm;
  int f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  int el2pel, inds_offset;
  int *inds;
  double *sqrt_tbl, *signs;
  int Ftm_stride, Ftm_offset;
  ssht_complex_double *Ftm;
  double *in_real;
  ssht_complex_double *out;
  double theta, ssign, elfactor;
  fftw_plan plan;
  double *thetas, *weights;
  double w;
  int spin = 0;

   // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
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

 // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing forward transform using GL sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_gl_forward_sov_real...");
  }

  // Compute weights and theta positions.
  thetas = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(thetas)
  weights = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(weights)
  ssht_sampling_gl_thetas_weights(thetas, weights, L);

  // Compute Fourier transform over phi, i.e. compute Ftm.
  Ftm = (ssht_complex_double*)calloc(L*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = L;
  Ftm_offset = 0;
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
      Ftm[t*Ftm_stride + m + Ftm_offset] =
	out[m] * 2.0 * SSHT_PI / (2.0*L-1.0);
  }
  free(in_real);
  free(out);
  fftw_destroy_plan(plan);

  // Compute flm.
  dlm1p1_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
  dl_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  inds_offset = 0;
  for (el=0; el<=L-1; el++) {
    for (m=0; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (t=0; t<=L-1; t++) {
    theta = thetas[t];
    w = weights[t];
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
      el2pel = el *el + el;
      for (m=0; m<=el; m++)
	inds[m + inds_offset] = el2pel + m;

      // Compute half dl line for given spin.
      ssht_dl_beta_kostelec_halfline_table(dlm1p1_line, dl_line,
					   theta, L, -spin, el,
					   sqrt_tbl, signs);

      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign
	  * elfactor
	  * w
	  * dl_line[m]
	  * Ftm[t*Ftm_stride + m + Ftm_offset];
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
  free(dlm1p1_line);
  free(dl_line);
  free(thetas);
  free(weights);
  free(Ftm);
  free(signs);
  free(sqrt_tbl);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Forward transform computed!");

}


//============================================================================
// DH algorithms
//============================================================================


/*!
 * Compute inverse transform using direct method with separation of
 * variables for DH sampling.
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_dh_inverse_sov(ssht_complex_double *f, const ssht_complex_double *flm,
			      int L, int spin, int verbosity) {

  int t, p, m, el, ind;
  int ftm_stride, ftm_offset, f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  double *sqrt_tbl, *signs;
  ssht_complex_double *ftm, *inout;
  double theta, ssign, elfactor;
  fftw_plan plan;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)

  // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  for (m=0; m<=L-1; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }
  ssign = signs[abs(spin)];

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing inverse transform using DH sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_dh_inverse_sov...");
  }

  // Compute ftm.
  ftm = (ssht_complex_double*)calloc((2*L)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(ftm)
  ftm_stride = 2*L-1;
  ftm_offset = L-1;
  dlm1p1_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
  dl_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  for (t=0; t<=2*L-1; t++) {
    theta = ssht_sampling_dh_t2theta(t, L);
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));

      // Compute dl line for given spin.
      ssht_dl_beta_kostelec_line_table(dlm1p1_line, dl_line,
      				       theta, L, -spin, el,
      				       sqrt_tbl, signs);
      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

      for (m=-el; m<=el; m++) {
	ssht_sampling_elm2ind(&ind, el, m);
	ftm[t*ftm_stride + m + ftm_offset] +=
	  ssign
	  * elfactor
	  * dl_line[m + L-1]
	  * flm[ind];
      }
    }
  }

  // Free dl memory.
  free(dlm1p1_line);
  free(dl_line);

  // Compute f.
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  f_stride = 2*L-1;
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
  for (t=0; t<=2*L-1; t++) {
    for (m=0; m<=L-1; m++)
      inout[m] = ftm[t*ftm_stride + m + ftm_offset];
    for (m=-(L-1); m<=-1; m++)
      inout[m+2*L-1] = ftm[t*ftm_stride + m + ftm_offset];
    fftw_execute_dft(plan, inout, inout);
    for (p=0; p<=2*L-2; p++)
      f[t*f_stride + p] = inout[p];
  }
  fftw_destroy_plan(plan);

  // Free memory.
  free(ftm);
  free(inout);
  free(signs);
  free(sqrt_tbl);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

}


/*!
 * Compute inverse transform of real scalar signal using direct method
 * with separation of variables for DH sampling (symmetries for real
 * signals are exploited).
 *
 * \param[out] f Function on sphere.
 * \param[in] flm Harmonic coefficients.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_dh_inverse_sov_real(double *f, const ssht_complex_double *flm,
				   int L, int verbosity) {

  int t, p, m, el, ind;
  int ftm_stride, ftm_offset, f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  double *sqrt_tbl, *signs;
  ssht_complex_double *ftm;
  ssht_complex_double *in;
  double *out;
  double theta, ssign, elfactor;
  fftw_plan plan;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)

  // Perform precomputations.
  for (el=0; el<=2*(L-1)+1; el++)
    sqrt_tbl[el] = sqrt((double)el);
  for (m=0; m<=L-1; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }
  ssign = signs[abs(spin)];

  // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing inverse transform using DH sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_dh_inverse_sov_real...");
  }

  // Compute ftm.
  ftm = (ssht_complex_double*)calloc(2*L*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(ftm)
  ftm_stride = L;
  ftm_offset = 0;
  dlm1p1_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
  dl_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  for (t=0; t<=2*L-1; t++) {
    theta = ssht_sampling_dh_t2theta(t, L);
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));

      // Compute half dl line for given spin.
      ssht_dl_beta_kostelec_halfline_table(dlm1p1_line, dl_line,
					   theta, L, -spin, el,
					   sqrt_tbl, signs);
      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

      for (m=0; m<=el; m++) {
	ssht_sampling_elm2ind(&ind, el, m);
	ftm[t*ftm_stride + m + ftm_offset] +=
	  ssign
	  * elfactor
	  * dl_line[m]
	  * flm[ind];
      }
    }
  }

  // Free dl memory.
  free(dlm1p1_line);
  free(dl_line);

  // Compute f.
  in = (ssht_complex_double*)calloc(L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in)
  out = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out)
  plan = fftw_plan_dft_c2r_1d(2*L-1, in, out, FFTW_MEASURE);
  f_stride = 2*L-1;
  for (t=0; t<=2*L-1; t++) {
    memcpy(in, &ftm[t*ftm_stride], L*sizeof(ssht_complex_double));
    fftw_execute_dft_c2r(plan, in, out);
    for (p=0; p<=2*L-2; p++)
      f[t*f_stride + p] = out[p];
  }
  fftw_destroy_plan(plan);

  // Free memory.
  free(ftm);
  free(in);
  free(out);
  free(signs);
  free(sqrt_tbl);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Inverse transform computed!");

}


/*!
 * Compute forward transform using Driscoll and Healy quadrature with
 * separation of variables.
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_dh_forward_sov(ssht_complex_double *flm, const ssht_complex_double *f,
			      int L, int spin, int verbosity) {

  int t, m, el, ind;
  int f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  int el2pel, inds_offset;
  int *inds;
  double *sqrt_tbl, *signs;
  int Ftm_stride, Ftm_offset;
  ssht_complex_double *Ftm, *inout;
  double theta, ssign, elfactor;
  fftw_plan plan;
  double w;

   // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
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

 // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing forward transform using DH sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_dh_forward_sov...");
  }

  // Compute Fourier transform over phi, i.e. compute Ftm.
  Ftm = (ssht_complex_double*)calloc((2*L)*(2*L-1), sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = 2*L-1;
  Ftm_offset = L-1;
  f_stride = 2*L-1;
  inout = (ssht_complex_double*)calloc(2*L-1, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (t=0; t<=2*L-1; t++) {
    memcpy(inout, &f[t*f_stride], f_stride*sizeof(ssht_complex_double));
    fftw_execute_dft(plan, inout, inout);
    for(m=0; m<=L-1; m++)
      Ftm[t*Ftm_stride + m + Ftm_offset] =
	inout[m] * 2.0 * SSHT_PI / (2.0*L-1.0) ;
    for(m=-(L-1); m<=-1; m++)
      Ftm[t*Ftm_stride + m + Ftm_offset] =
	inout[m+2*L-1] * 2.0 * SSHT_PI / (2.0*L-1.0);
  }
  fftw_destroy_plan(plan);

  // Compute flm.
  dlm1p1_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
  dl_line = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  inds_offset = L-1;
  for (el=0; el<=L-1; el++) {
    for (m=-el; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (t=0; t<=2*L-1; t++) {
    theta = ssht_sampling_dh_t2theta(t, L);
    w = ssht_sampling_weight_dh(theta, L);
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
      el2pel = el *el + el;
      for (m=-el; m<=el; m++)
	inds[m + inds_offset] = el2pel + m;

      // Compute dl line for given spin.
      ssht_dl_beta_kostelec_line_table(dlm1p1_line, dl_line,
      				       theta, L, -spin, el,
      				       sqrt_tbl, signs);
      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

      for (m=-el; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign
	  * elfactor
	  * w
	  * dl_line[m+L-1]
	  * Ftm[t*Ftm_stride + m + Ftm_offset];
      }
    }
  }

  // Free memory.
  free(dlm1p1_line);
  free(dl_line);
  free(Ftm);
  free(inout);
  free(signs);
  free(sqrt_tbl);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Forward transform computed!");

}


/*!
 * Compute forward transform of real scalar signal using Driscoll and
 * Healy quadrature with separation of variables (symmetries for real
 * signals are exploited).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosity flag in range [0,5].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_core_dh_forward_sov_real(ssht_complex_double *flm, const double *f,
				   int L, int verbosity) {

  int t, m, el, ind, ind_nm;
  int f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  int el2pel, inds_offset;
  int *inds;
  double *sqrt_tbl, *signs;
  int Ftm_stride, Ftm_offset;
  ssht_complex_double *Ftm;
  double *in_real;
  ssht_complex_double *out;
  double theta, ssign, elfactor;
  fftw_plan plan;
  double w;
  int spin = 0;

   // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
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

 // Print messages depending on verbosity level.
  if (verbosity > 0) {
    printf("%s%s\n", SSHT_PROMPT,
	   "Computing forward transform using GL sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (",
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s%s\n", SSHT_PROMPT,
	     "Using routine ssht_core_gl_forward_sov_real...");
  }

  // Compute Fourier transform over phi, i.e. compute Ftm.
  Ftm = (ssht_complex_double*)calloc(2*L*L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = L;
  Ftm_offset = 0;
  f_stride = 2*L-1;
  in_real = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in_real)
  out = (ssht_complex_double*)calloc(L, sizeof(ssht_complex_double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out)
  plan = fftw_plan_dft_r2c_1d(2*L-1, in_real, out, FFTW_MEASURE);
  for (t=0; t<=2*L-1; t++) {
    memcpy(in_real, &f[t*f_stride], f_stride*sizeof(double));
    fftw_execute_dft_r2c(plan, in_real, out);
    for(m=0; m<=L-1; m++)
      Ftm[t*Ftm_stride + m + Ftm_offset] =
	out[m] * 2.0 * SSHT_PI / (2.0*L-1.0);
  }
  free(in_real);
  free(out);
  fftw_destroy_plan(plan);

  // Compute flm.
  dlm1p1_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
  dl_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  inds_offset = 0;
  for (el=0; el<=L-1; el++) {
    for (m=0; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (t=0; t<=2*L-1; t++) {
    theta = ssht_sampling_dh_t2theta(t, L);
    w = ssht_sampling_weight_dh(theta, L);
    for (el=abs(spin); el<=L-1; el++) {
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
      el2pel = el *el + el;
      for (m=0; m<=el; m++)
	inds[m + inds_offset] = el2pel + m;

      // Compute half dl line for given spin.
      ssht_dl_beta_kostelec_halfline_table(dlm1p1_line, dl_line,
					   theta, L, -spin, el,
					   sqrt_tbl, signs);

      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign
	  * elfactor
	  * w
	  * dl_line[m]
	  * Ftm[t*Ftm_stride + m + Ftm_offset];
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
  free(dlm1p1_line);
  free(dl_line);
  free(Ftm);
  free(signs);
  free(sqrt_tbl);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0)
    printf("%s%s", SSHT_PROMPT, "Forward transform computed!");

}
