



#include <stdio.h>
#include <stdlib.h>
#include <string.h> // for memcpy
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "ssht_types.h"
#include "ssht_error.h"
#include "ssht_dl.h"
#include "ssht_sampling.h"
#include "ssht_core.h"




void ssht_core_mw_inverse_sov_sym(complex double *f, complex double *flm, 
				  int L, int spin, int verbosity) {

  int el, m, mm, ind, t, p;
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
	   "Computing inverse transform using McEwen and Wiaux sampling with ");
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
  dl_offset = ssht_dl_get_mmoffset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_mmstride(L, SSHT_DL_QUARTER);   
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
      ssht_dl_halfpi_trapani_quarter_table(dl, L,
    					  SSHT_DL_QUARTER,
    					  el, sqrt_tbl);
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










// note flm must be initialised with zeros.
void ssht_core_mw_forward_sov_conv_sym(complex double *flm, complex double *f, 
				       int L, int spin, int verbosity) {

  int el, m, mm, ind, t, r;
  int eltmp;
  double *sqrt_tbl, *signs;
  int el2pel, inds_offset;
  int *inds;
  double ssign, elfactor;
  fftw_plan plan, plan_bwd, plan_fwd;
  //complex double *in, *out;
  complex double *inout, *in_bwd, *out_bwd, *in_fwd, *out_fwd;
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
	   "Computing forward transform using McEwen and Wiaux sampling with ");
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
//**TODO: use memcpy.
  //memcpy(&wr[w_offset+1], &w[0], (2*L-2)*sizeof(complex double));
  //memcpy(&wr[0], &w[w_offset], (2*L-1)*sizeof(complex double));
  fftw_execute_dft(plan_bwd, inout, inout);
  for (mm=0; mm<=2*L-2; mm++) 
    wr[mm + w_offset] = inout[mm - 2*(L-1) + w_offset];
  for (mm=-2*(L-1); mm<=-1; mm++) 
    wr[mm + w_offset] = inout[mm + 2*(L-1) + 1 + w_offset];
//**TODO: use memcpy.
  //memcpy(&wr[w_offset], &w[0], (2*L-1)*sizeof(complex double));
  //memcpy(&wr[0], &w[w_offset+1], (2*L-2)*sizeof(complex double));

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
//      Gmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] = 
//	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);
      Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] = 
	Fmm_pad[mm + w_offset] * 2.0 * SSHT_PI / (4.0*L-3.0);


  }
  fftw_destroy_plan(plan_bwd);
  fftw_destroy_plan(plan_fwd);

  // Compute flm.
  dl = ssht_dl_calloc(L, SSHT_DL_QUARTER);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  dl_offset = ssht_dl_get_mmoffset(L, SSHT_DL_QUARTER);
  dl_stride = ssht_dl_get_mmstride(L, SSHT_DL_QUARTER); 
  for (el=0; el<=L-1; el++) {
    for (m=-el; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      flm[ind] = 0.0;
    }
  }
  for (el=abs(spin); el<=L-1; el++) {

    // Compute Wigner plane.
    /* if (el!=0 && el==abs(spin)) { */
    /*   for(eltmp=0; eltmp<=abs(spin); eltmp++) { */
    /* 	ssht_dl_halfpi_trapani_eighth_table(dl, L, */
    /* 					    SSHT_DL_HALF, */
    /* 					    eltmp, sqrt_tbl); */
    /*   } */
    /*   ssht_dl_halfpi_trapani_fill_eighth2righthalf_table(dl, L, */
    /* 							 SSHT_DL_HALF, */
    /* 							 el, signs); */
    /* } */
    /* else { */
    /*   ssht_dl_halfpi_trapani_eighth_table(dl, L, */
    /* 					  SSHT_DL_HALF, */
    /* 					  el, sqrt_tbl); */
    /*   ssht_dl_halfpi_trapani_fill_eighth2righthalf_table(dl, L, */
    /* 							 SSHT_DL_HALF, */
    /* 							 el, signs); */
    /* } */

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
      ssht_dl_halfpi_trapani_quarter_table(dl, L,
    					  SSHT_DL_QUARTER,
    					  el, sqrt_tbl);
    }

    // Compute flm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    /* for (m=-el; m<=el; m++) { */
    /*   ssht_sampling_elm2ind(&ind, el, m); */

    /*   /\* flm[ind] += *\/ */
    /*   /\* 	ssign  *\/ */
    /*   /\* 	* elfactor *\/ */
    /*   /\* 	* cexp(I*SSHT_PION2*(m+spin)) *\/ */
    /*   /\* 	* dl[0*dl_stride + m + dl_offset] *\/ */
    /*   /\* 	* dl[0*dl_stride - spin + dl_offset] *\/ */
    /*   /\* 	* Gmm[(m+Fmm_offset)*Fmm_stride + 0 + Fmm_offset]; *\/ */

    /*   flm[ind] += */
    /* 	ssign  */
    /* 	* elfactor */
    /* 	* cexp(I*SSHT_PION2*(m+spin)) */
    /* 	* signs[el] * dl[0*dl_stride + abs(m) + dl_offset] */
    /* 	* dl[0*dl_stride - spin + dl_offset] */
    /* 	* Gmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset]; */

    /*   /\* for (mm=1; mm<=el; mm++) *\/ */
    /*   /\* 	flm[ind] += *\/ */
    /*   /\* 	  ssign  *\/ */
    /*   /\* 	  * elfactor *\/ */
    /*   /\* 	  * cexp(I*SSHT_PION2*(m+spin)) *\/ */
    /*   /\* 	  * dl[mm*dl_stride + m + dl_offset] *\/ */
    /*   /\* 	  * dl[mm*dl_stride - spin + dl_offset] *\/ */
    /*   /\* 	  * ( Gmm[(m+Fmm_offset)*Fmm_stride + mm + Fmm_offset] *\/ */
    /*   /\* 	      + signs[abs(m)] * ssign *\/ */
    /*   /\* 	      * Gmm[(m+Fmm_offset)*Fmm_stride - mm + Fmm_offset]); *\/ */

    /*   for (mm=1; mm<=el; mm++) { */
    /* 	elmmsign = m < 0 ? signs[el] * signs[mm] : 1.0; */
    /* 	elssign = spin <= 0 ? 1.0 : elmmsign; */

    /* 	flm[ind] += */
    /* 	  ssign  */
    /* 	  * elfactor */
    /* 	  * cexp(I*SSHT_PION2*(m+spin)) */
    /* 	  * elmmsign * dl[mm*dl_stride + abs(m) + dl_offset] */
    /* 	  * dl[mm*dl_stride - spin + dl_offset] */
    /* 	  * ( Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] */
    /* 	      + signs[abs(m)] * ssign */
    /* 	      * Gmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]); */
    /*   } */
    /* }  */ 


    el2pel = el *el + el;    
    for (m=-el; m<=el; m++)
      inds[m + inds_offset] = el2pel + m; 
    for (m=-el; m<=el; m++) {
      // mm = 0
      ind = inds[m + inds_offset];
      flm[ind] +=
	ssign 
	* elfactor
	* expsm[m + exps_offset]
	* signs[el] * dl[0*dl_stride + abs(m) + dl_offset]
	* dl[0*dl_stride - spin + dl_offset]
	* Gmm[(0+Fmm_offset)*Fmm_stride + m + Fmm_offset];
    }

    for (mm=1; mm<=el; mm++) {
      elmmsign = m < 0 ? signs[el] * signs[mm] : 1.0;
      elssign = spin <= 0 ? 1.0 : signs[el] * signs[mm];

      for (m=-el; m<=el; m++) {
	ind = inds[m + inds_offset];
	flm[ind] +=
	  ssign 
	  * elfactor
	  * expsm[m + exps_offset]
	  * elmmsign * dl[mm*dl_stride + abs(m) + dl_offset]
	  * elssign * dl[mm*dl_stride - spinneg + dl_offset]
	  * ( Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]
	      + signs[abs(m)] * ssign
	      * Gmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]);
      }

      /* for (m=-el; m<=-1; m++) { */
      /*   ind = inds[m + inds_offset]; */
      /*   flm[ind] += */
      /*     ssign  */
      /*     * elfactor */
      /*     * expsm[m + exps_offset] */
      /*     * elmmsign * dl[mm*dl_stride - m + dl_offset] */
      /*     * elssign * dl[mm*dl_stride - spinneg + dl_offset] */
      /*     * ( Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] */
      /* 	+ signs[-m] * ssign */
      /* 	* Gmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]); */
      /* } */
      /* for (m=0; m<=el; m++) { */
      /*   ind = inds[m + inds_offset]; */
      /*   flm[ind] += */
      /*     ssign  */
      /*     * elfactor */
      /*     * expsm[m + exps_offset] */
      /*     * elmmsign * dl[mm*dl_stride + m + dl_offset] */
      /*     * elssign * dl[mm*dl_stride - spinneg + dl_offset] */
      /*     * ( Gmm[(mm+Fmm_offset)*Fmm_stride + m + Fmm_offset] */
      /* 	+ signs[m] * ssign */
      /* 	* Gmm[(-mm+Fmm_offset)*Fmm_stride + m + Fmm_offset]); */
      /* } */
    }  


  }

  // Free dl memory.
  free(dl);


  /* for (el=abs(spin); el<=L-1; el++) { */
  /*   for (m=-el; m<=el; m++) { */
  /*     ssht_sampling_elm2ind(&ind, el, m); */
  /*     flm[ind] = 1.0 + I * 2.0; */
  /*   } */
  /* } */




  
  free(Fmt);
  free(Fmm);
  //  free(out);
  free(inout);
  //  free(inout_bwd);

  free(w);
  free(wr);
  free(Fmm_pad);
  free(tmp_pad);
  free(Gmm);
  

  free(expsm);
  free(expsmm);
  free(inds);

  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s %s", SSHT_PROMPT, "Forward transform computed!");  

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs); 

}
