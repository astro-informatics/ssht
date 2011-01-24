



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#include "ssht_types.h"
#include "ssht_error.h"
#include "ssht_dl.h"
#include "ssht_core.h"




void ssht_core_mw_inverse_sov_sym(complex double *f, complex double *flm, 
				  int L, int spin, int verbosity) {


  int el, m, mm, ind, t, p;
  int eltmp;
  double *sqrt_tbl, *signs;
  double ssign, elfactor;
  double *dl;
  int dl_offset, dl_stride;
  complex double *Fmm, *fext;
  int Fmm_offset, Fmm_stride, fext_stride;
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
  dl = ssht_dl_calloc(L, SSHT_DL_HALF);
  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  dl_offset = ssht_dl_get_mmoffset(L, SSHT_DL_HALF);
  dl_stride = ssht_dl_get_mmstride(L, SSHT_DL_HALF);   
  for (el=abs(spin); el<=L-1; el++) {

    // Compute Wigner plane.
    if (el!=0 && el==abs(spin)) {
      for(eltmp=0; eltmp<=abs(spin); eltmp++) {
	ssht_dl_halfpi_trapani_eighth_table(dl, L, 
					    SSHT_DL_HALF,
					    eltmp, sqrt_tbl);
      }
      ssht_dl_halfpi_trapani_fill_eighth2righthalf_table(dl, L,
							 SSHT_DL_HALF,
							 el, signs);
    }
    else {
      ssht_dl_halfpi_trapani_eighth_table(dl, L, 
					  SSHT_DL_HALF,
					  el, sqrt_tbl);
      ssht_dl_halfpi_trapani_fill_eighth2righthalf_table(dl, L,
							SSHT_DL_HALF,
							el, signs);
    }

    // Compute Fmm.
    elfactor = sqrt((double)(2.0*el+1.0)/(4.0*SSHT_PI));
    for (m=-el; m<=el; m++) {
      ssht_sampling_elm2ind(&ind, el, m);
      for (mm=0; mm<=el; mm++) {
	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] = 
	  Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset]
	  + ssign 
	  * elfactor
	  * cexp(-I*SSHT_PION2*(m+spin))
	  * dl[mm*dl_stride + m + dl_offset] 
	  * dl[mm*dl_stride - spin + dl_offset]
	  * flm[ind];
      }
    }  

  }

  // Free dl memory.
  free(dl);

  // Use symmetry to compute Fmm for negative mm.
  for (m=-(L-1); m<=L-1; m++) {
    for (mm=-(L-1); mm<=-1; mm++) {
      Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] = 
	signs[abs(m)] * ssign 
	* Fmm[(m + Fmm_offset)*Fmm_stride - mm + Fmm_offset];
    }
  }

  // Apply phase modulation to account for sampling offset.
  for (m=-(L-1); m<=L-1; m++) {
    for (mm=-(L-1); mm<=L-1; mm++) {
      Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset] = 
	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset]
	* cexp(I*mm*SSHT_PI/(2.0*L-1.0));
    }
  }

  // Allocate space for function values.
  fext = (complex double*)calloc((2*L-1)*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(fext)
  fext_stride = 2*L-1;

  // Plan fftw before intialise memory.
  plan = fftw_plan_dft_2d(2*L-1, 2*L-1, fext, fext, 
			       FFTW_BACKWARD, FFTW_ESTIMATE);

  // Apply spatial shift.
  for (m=0; m<=L-1; m++)
    for (mm=0; mm<=L-1; mm++)
      fext[m*fext_stride + mm] = 
	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset];
  for (m=-(L-1); m<=-1; m++)
    for (mm=0; mm<=L-1; mm++)
      fext[(m+2*L-1)*fext_stride + mm] = 
	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset];
  for (m=0; m<=L-1; m++)
    for (mm=-(L-1); mm<=-1; mm++)
      fext[m*fext_stride + mm + 2*L-1] = 
	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset];
  for (m=-(L-1); m<=-1; m++)
    for (mm=-(L-1); mm<=-1; mm++)
      fext[(m+2*L-1)*fext_stride + mm + 2*L-1] = 
	Fmm[(m + Fmm_offset)*Fmm_stride + mm + Fmm_offset];

  // Free Fmm memory.
  free(Fmm);

  // Perform 2D FFT.  
  fftw_execute(plan);
  fftw_destroy_plan(plan);
  
  // Extract f from version of f extended to the torus (fext).
  for (t=0; t<=L-1; t++)
    for (p=0; p<=2*L-2; p++)
      f[t*fext_stride + p] = fext[t*fext_stride + p];
 
  // Free fext memory.
  free(fext);
  
  // Print finished if verbosity set.
  if (verbosity > 0) 
    printf("%s %s", SSHT_PROMPT, "Inverse transform computed!");  

  // Free precomputation memory.
  free(sqrt_tbl);
  free(signs); 

}
