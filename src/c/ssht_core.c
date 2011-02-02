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
#include <complex.h>  // Must be before fftw3.h
#include <fftw3.h>

#include "ssht_types.h"
#include "ssht_error.h"
#include "ssht_dl.h"
#include "ssht_sampling.h"
#include "ssht_core.h"


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
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_mw_inverse_sov_sym(complex double *f, complex double *flm, 
				  int L, int spin, int verbosity) {

  int el, m, mm, ind;
  //int t, p;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  complex double mmfactor;
  double *dl;
  int dl_offset, dl_stride;
  complex double *exps;
  int exps_offset;
  double elmmsign, elssign;
  int spinneg;
  complex double *Fmm, *fext;
  int Fmm_offset, Fmm_stride, fext_stride;
  fftw_plan plan;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  exps = (complex double*)calloc(2*L-1, sizeof(complex double));
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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_mw_inverse_sov_sym...");
  }

  // Compute Fmm.
  Fmm = (complex double*)calloc((2*L-1)*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_offset = L-1;
  Fmm_stride = 2*L-1;    
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);   
  inds_offset = L-1;
  for (el=abs(spin); el<=L-1; el++) {

    // Compute Wigner plane.
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
  fext = (complex double*)calloc((2*L-1)*(2*L-1), sizeof(complex double));
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
  memcpy(f, fext, L*(2*L-1)*sizeof(complex double));
  /* Memcpy equivalent to:
  for (t=0; t<=L-1; t++)
    for (p=0; p<=2*L-2; p++)
      f[t*fext_stride + p] = fext[t*fext_stride + p];
  */

  // Free fext memory.
  free(fext);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s %s", SSHT_PROMPT, "Inverse transform computed!");  

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs); 
  free(exps);
  free(inds);

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
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_mw_inverse_sov_sym_real(double *f, complex double *flm, 
				       int L, int verbosity) {

  int el, m, mm, ind;
  //int t, p;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  complex double mmfactor;
  double *dl;
  int dl_offset, dl_stride;
  complex double *exps;
  int exps_offset;
  double elmmsign, elssign;
  int spinneg;
  complex double *Fmm, *Fmm_shift;
  double *fext_real;
  int Fmm_offset, Fmm_stride, fext_stride;
  fftw_plan plan;
  int spin = 0;

  // Allocate memory.
  sqrt_tbl = (double*)calloc(2*(L-1)+2, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(sqrt_tbl)
  signs = (double*)calloc(L+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(signs)
  exps = (complex double*)calloc(2*L-1, sizeof(complex double));
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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_mw_inverse_sov_sym_real...");
  }

  // Compute Fmm.
  Fmm = (complex double*)calloc((2*L-1)*L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_offset = L-1;
  Fmm_stride = L;    
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  dl_offset = ssht_dl_get_offset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_stride(L, SSHT_DL_QUARTER);   
  inds_offset = L-1;
  for (el=abs(spin); el<=L-1; el++) {

    // Compute Wigner plane.
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
  Fmm_shift = (complex double*)calloc((2*L-1)*L, sizeof(complex double));
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
  fext_stride = 2*L-1;

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
    printf("%s %s", SSHT_PROMPT, "Inverse transform computed!");  

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs); 
  free(exps);
  free(inds);

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
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_mwdirect_inverse(complex double *f, complex double *flm, 
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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
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
    printf("%s %s", SSHT_PROMPT, "Inverse transform computed!");  

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
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_mwdirect_inverse_sov(complex double *f, complex double *flm, 
				     int L, int spin, int verbosity) {

  int t, p, m, el, ind;
  int ftm_stride, ftm_offset, f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  double *sqrt_tbl, *signs;
  complex double *ftm, *inout;
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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing inverse transform using MW sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_mwdirect_inverse_sov...");
  }

  // Compute ftm.
  ftm = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
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
  inout = (complex double*)calloc(2*L-1, sizeof(complex double));
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
    printf("%s %s", SSHT_PROMPT, "Inverse transform computed!");  

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
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_mw_forward_sov_conv_sym(complex double *flm, complex double *f, 
				       int L, int spin, int verbosity) {

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
  Fmm = (complex double*)calloc((2*L-1)*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Fmm)
  Fmm_stride = 2*L-1;
  Fmm_offset = L-1;
  for (m=-(L-1); m<=L-1; m++) {
    memcpy(inout, &Fmt[(m+Fmt_offset)*Fmt_stride], Fmt_stride*sizeof(complex double));
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

  // Compute flm.
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
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
	* Gmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
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
	* Gmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
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
      	  * ( Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]
      	      + signs[-m] * ssign
      	      * Gmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]);
      }
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign 
	  * elfactor
	  * expsm[m + exps_offset]
	  * dl[mm*dl_stride + m + dl_offset]
	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
	  * ( Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]
	      + signs[m] * ssign
	      * Gmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]);
      }

    }  

  }

  // Free memory.
  free(dl);
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
 * Compute forward transform for MW method of real scalar signal using
 * separation of variables, fast Fourier transforms, performing
 * convolution with weights as product in transformed space and
 * exploiting all symmetries (including additional symmetries for real
 * signals).
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_mw_forward_sov_conv_sym_real(complex double *flm, double *f, 
					    int L, int verbosity) {

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
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_gl_inverse_sov(complex double *f, complex double *flm, 
			      int L, int spin, int verbosity) {

  int t, p, m, el, ind;
  int ftm_stride, ftm_offset, f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  double *sqrt_tbl, *signs;
  complex double *ftm, *inout;
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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing inverse transform using GL sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_gl_inverse_sov...");
  }

  // Compute weights and theta positions.
  thetas = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(thetas)
  weights = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(weights)
  ssht_sampling_gl_thetas_weights(thetas, weights, L);

  // Compute ftm.
  ftm = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
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
  inout = (complex double*)calloc(2*L-1, sizeof(complex double));
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
    printf("%s %s", SSHT_PROMPT, "Inverse transform computed!");  

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
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_gl_inverse_sov_real(double *f, complex double *flm, 
				   int L, int verbosity) {

  int t, p, m, el, ind;
  int ftm_stride, ftm_offset, f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  double *sqrt_tbl, *signs;
  complex double *ftm;
//, *inout;
  double theta, ssign, elfactor;
  fftw_plan plan;
  double *thetas, *weights;
  int spin = 0;


complex double *in;
double *out;


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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing inverse transform using GL sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_gl_inverse_sov_real...");
  }

  // Compute weights and theta positions.
  thetas = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(thetas)
  weights = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(weights)
  ssht_sampling_gl_thetas_weights(thetas, weights, L);

  // Compute ftm.
//  ftm = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
  ftm = (complex double*)calloc(L*L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(ftm)  
//  ftm_stride = 2*L-1;
  ftm_stride = L;
//  ftm_offset = L-1;
  ftm_offset = 0;
//  dlm1p1_line = (double*)calloc(2*L-1, sizeof(double));
  dlm1p1_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
//  dl_line = (double*)calloc(2*L-1, sizeof(double));
  dl_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  for (t=0; t<=L-1; t++) {
    theta = thetas[t];
    for (el=abs(spin); el<=L-1; el++) {	
      elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));

      // Compute dl line for given spin.
      /* ssht_dl_beta_kostelec_line_table(dlm1p1_line, dl_line, */
      /* 				       theta, L, -spin, el, */
      /* 				       sqrt_tbl, signs); */

      ssht_dl_beta_kostelec_halfline_table(dlm1p1_line, dl_line,
					   theta, L, -spin, el,
					   sqrt_tbl, signs);
      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;
    
//      for (m=-el; m<=el; m++) {	
      for (m=0; m<=el; m++) {	
	ssht_sampling_elm2ind(&ind, el, m);
	ftm[t*ftm_stride + m + ftm_offset] +=
	  ssign 
	  * elfactor
//	  * dl_line[m + L-1]
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

//  inout = (complex double*)calloc(2*L-1, sizeof(complex double));
//  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
//  f_stride = 2*L-1;
//  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_BACKWARD, FFTW_MEASURE);
//  for (t=0; t<=L-1; t++) {
//    for (m=0; m<=L-1; m++)
//      inout[m] = ftm[t*ftm_stride + m + ftm_offset];
//    for (m=-(L-1); m<=-1; m++)
//      inout[m+2*L-1] = ftm[t*ftm_stride + m + ftm_offset];
//    fftw_execute_dft(plan, inout, inout);
//    for (p=0; p<=2*L-2; p++)
//      f[t*f_stride + p] = inout[p];
//  }

  in = (complex double*)calloc(L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(in)
  out = (double*)calloc(2*L-1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(out)
  plan = fftw_plan_dft_c2r_1d(2*L-1, in, out, FFTW_MEASURE);
  f_stride = 2*L-1;
  for (t=0; t<=L-1; t++) {
    memcpy(in, &ftm[t*ftm_stride], L*sizeof(complex double));
    fftw_execute_dft_c2r(plan, in, out);
    for (p=0; p<=2*L-2; p++)
      f[t*f_stride + p] = out[p];
  }
  fftw_destroy_plan(plan);

  // Free memory.  
  free(ftm);
//  free(inout);
  free(in);
  free(out);
  free(signs);
  free(sqrt_tbl);
  
  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s %s", SSHT_PROMPT, "Inverse transform computed!");  

}


/*!  
 * Compute forward transform using Gauss-Legendgre quadrature with separation of
 * variables.
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_gl_forward_sov(complex double *flm, complex double *f, 
			      int L, int spin, int verbosity) {

  int t, m, el, ind;
  int f_stride;  
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  int el2pel, inds_offset;
  int *inds;
  double *sqrt_tbl, *signs;
  int Ftm_stride, Ftm_offset;
  complex double *Ftm, *inout;
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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing forward transform using GL sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_gl_forward_sov...");
  }

  // Compute weights and theta positions.
  thetas = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(thetas)
  weights = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(weights)
  ssht_sampling_gl_thetas_weights(thetas, weights, L);

  // Compute Fourier transform over phi, i.e. compute Ftm.
  Ftm = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = 2*L-1;
  Ftm_offset = L-1;
  f_stride = 2*L-1;
  inout = (complex double*)calloc(2*L-1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (t=0; t<=L-1; t++) {
    memcpy(inout, &f[t*f_stride], f_stride*sizeof(double complex));
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
    printf("%s %s", SSHT_PROMPT, "Forward transform computed!");  

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
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_gl_forward_sov_real(complex double *flm, double *f, 
				   int L, int verbosity) {

  int t, m, el, ind, ind_nm;
  int f_stride;  
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  int el2pel, inds_offset;
  int *inds;
  double *sqrt_tbl, *signs;
  int Ftm_stride, Ftm_offset;
  //  complex double *Ftm, *inout;

  complex double *Ftm;
double *in_real;
complex double *out;

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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing forward transform using GL sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", TRUE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_gl_forward_sov_real...");
  }

  // Compute weights and theta positions.
  thetas = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(thetas)
  weights = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(weights)
  ssht_sampling_gl_thetas_weights(thetas, weights, L);

  // Compute Fourier transform over phi, i.e. compute Ftm.
  /* Ftm = (complex double*)calloc(L*(2*L-1), sizeof(complex double)); */
  /* SSHT_ERROR_MEM_ALLOC_CHECK(Ftm) */
  /* Ftm_stride = 2*L-1; */
  /* Ftm_offset = L-1; */
  /* f_stride = 2*L-1; */
  /* inout = (complex double*)calloc(2*L-1, sizeof(complex double)); */
  /* SSHT_ERROR_MEM_ALLOC_CHECK(inout) */
  /* plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE); */
  /* for (t=0; t<=L-1; t++) { */
  /*   memcpy(inout, &f[t*f_stride], f_stride*sizeof(double complex)); */
  /*   fftw_execute_dft(plan, inout, inout); */
  /*   for(m=0; m<=L-1; m++)  */
  /*     Ftm[t*Ftm_stride + m + Ftm_offset] =  */
  /* 	inout[m] * 2.0 * SSHT_PI / (2.0*L-1.0) ; */
  /*   for(m=-(L-1); m<=-1; m++)  */
  /*     Ftm[t*Ftm_stride + m + Ftm_offset] =  */
  /* 	inout[m+2*L-1] * 2.0 * SSHT_PI / (2.0*L-1.0); */
  /* } */
  /* fftw_destroy_plan(plan); */

  // Compute Fourier transform over phi, i.e. compute Ftm.
  Ftm = (complex double*)calloc(L*L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = L;
  Ftm_offset = 0;
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
      Ftm[t*Ftm_stride + m + Ftm_offset] = 
	out[m] * 2.0 * SSHT_PI / (2.0*L-1.0);
  }
  free(in_real);
  free(out);
  fftw_destroy_plan(plan);

  // Compute flm.
//  dlm1p1_line = (double*)calloc(2*L-1, sizeof(double));
  dlm1p1_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dlm1p1_line)
//  dl_line = (double*)calloc(2*L-1, sizeof(double));
  dl_line = (double*)calloc(L, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dl_line)
  inds_offset = L-1;
  for (el=0; el<=L-1; el++) {
//    for (m=-el; m<=el; m++) {
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
//      for (m=-el; m<=el; m++)
      for (m=0; m<=el; m++)
	inds[m + inds_offset] = el2pel + m; 

      // Compute dl line for given spin.
      /* ssht_dl_beta_kostelec_line_table(dlm1p1_line, dl_line, */
      /* 				       theta, L, -spin, el, */
      /* 				       sqrt_tbl, signs); */

      ssht_dl_beta_kostelec_halfline_table(dlm1p1_line, dl_line,
					   theta, L, -spin, el,
					   sqrt_tbl, signs);

      // Switch current and previous dls.
      dl_ptr = dl_line;
      dl_line = dlm1p1_line;
      dlm1p1_line = dl_ptr;

//      for (m=-el; m<=el; m++) {
      for (m=0; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] += 
	  ssign
	  * elfactor
	  * w
//	  * dl_line[m+L-1]
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


//  free(inout);

  free(signs);
  free(sqrt_tbl);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s %s", SSHT_PROMPT, "Forward transform computed!");  

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
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_dh_inverse_sov(complex double *f, complex double *flm, 
			      int L, int spin, int verbosity) {

  int t, p, m, el, ind;
  int ftm_stride, ftm_offset, f_stride;
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  double *sqrt_tbl, *signs;
  complex double *ftm, *inout;
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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing inverse transform using DH sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_dh_inverse_sov...");
  }

  // Compute ftm.
  ftm = (complex double*)calloc((2*L)*(2*L-1), sizeof(complex double));
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
  inout = (complex double*)calloc(2*L-1, sizeof(complex double));
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
    printf("%s %s", SSHT_PROMPT, "Inverse transform computed!");  

}


/*!  
 * Compute forward transform using Driscoll and Healy quadrature with
 * separation of variables.
 *
 * \param[out] flm Harmonic coefficients.
 * \param[in] f Function on sphere.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] verbosity Verbosiity flag in range [0,5].
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_core_dh_forward_sov(complex double *flm, complex double *f, 
			      int L, int spin, int verbosity) {

  int t, m, el, ind;
  int f_stride;  
  double *dlm1p1_line,  *dl_line;
  double *dl_ptr;
  int el2pel, inds_offset;
  int *inds;
  double *sqrt_tbl, *signs;
  int Ftm_stride, Ftm_offset;
  complex double *Ftm, *inout;
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
    printf("%s %s\n", SSHT_PROMPT, 
	   "Computing forward transform using DH sampling with ");
    printf("%s%s%d%s%d%s\n", SSHT_PROMPT, "parameters  (L,spin,reality) = (", 
	   L, ",", spin, ", FALSE)");
    if (verbosity > 1)
      printf("%s %s\n", SSHT_PROMPT, 
	     "Using routine ssht_core_dh_forward_sov...");
  }

  // Compute Fourier transform over phi, i.e. compute Ftm.
  Ftm = (complex double*)calloc((2*L)*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(Ftm)
  Ftm_stride = 2*L-1;
  Ftm_offset = L-1;
  f_stride = 2*L-1;
  inout = (complex double*)calloc(2*L-1, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(inout)
  plan = fftw_plan_dft_1d(2*L-1, inout, inout, FFTW_FORWARD, FFTW_MEASURE);
  for (t=0; t<=2*L-1; t++) {
    memcpy(inout, &f[t*f_stride], f_stride*sizeof(double complex));
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
    printf("%s %s", SSHT_PROMPT, "Forward transform computed!");  

}
