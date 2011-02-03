/*! 
 * \file ssht_test.c
 * Applies SSHT algorithms to perform inverse and forward spherical
 * harmonic transforms (respectively) to check that the original
 * signal is reconstructed exactly (to numerical precision).  Test is
 * performed on a random signal with harmonic coefficients uniformly
 * sampled from (-1,1).
 *
 * Usage: ssht_test B spin, e.g. ssht_test 64 2
 *
 * \author Jason McEwen
 */

#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <fftw3.h>  // Must be before fftw3.h
#include <time.h>

#include <ssht.h>

#define NREPEAT 3
#define MAX(a,b) ((a) > (b) ? (a) : (b))

double ran2_dp(int idum);
void ssht_test_gen_flm_complex(complex double *flm, int L, int spin, int seed);
void ssht_test_gen_flm_real(complex double *flm, int L, int seed);


int main(int argc, char *argv[]) {

  complex double *flm_orig, *flm_syn;
  complex double *f_mw, *f_gl, *f_dh;
  double *f_mw_real, *f_gl_real, *f_dh_real;


complex double *f_mw_ss, *f_mw_ss2;



  int L = 128;
  int spin = 0;
  int irepeat;
  int seed = 1;
  int verbosity = 0;
  int i;
  double tmp;

  clock_t time_start, time_end;
  double errors_mw[NREPEAT];
  double errors_gl[NREPEAT];
  double errors_dh[NREPEAT];
  double errors_mw_real[NREPEAT];
  double errors_gl_real[NREPEAT];
  double errors_dh_real[NREPEAT];
  double durations_forward_mw[NREPEAT];
  double durations_inverse_mw[NREPEAT];
  double durations_forward_gl[NREPEAT];
  double durations_inverse_gl[NREPEAT];
  double durations_forward_dh[NREPEAT];
  double durations_inverse_dh[NREPEAT];
  double durations_forward_mw_real[NREPEAT];
  double durations_inverse_mw_real[NREPEAT];
  double durations_forward_gl_real[NREPEAT];
  double durations_inverse_gl_real[NREPEAT];
  double durations_forward_dh_real[NREPEAT];
  double durations_inverse_dh_real[NREPEAT];

  // Parse problem sizes.
  L = atoi(argv[1]);
  spin = atoi(argv[2]);

  // Allocate memory.
  flm_orig = (complex double*)calloc(L*L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(flm_orig)
  flm_syn = (complex double*)calloc(L*L, sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(flm_syn)
  f_mw = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_mw)
  f_gl = (complex double*)calloc(L*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_gl)
  f_dh = (complex double*)calloc((2*L)*(2*L-1), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_dh)
  f_mw_real = (double*)calloc(L*(2*L-1), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_mw_real)
  f_gl_real = (double*)calloc(L*(2*L-1), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_gl_real)
  f_dh_real = (double*)calloc((2*L)*(2*L-1), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(f_dh_real)



f_mw_ss = (complex double*)calloc((L+1)*(2*L-1), sizeof(complex double));
SSHT_ERROR_MEM_ALLOC_CHECK(f_mw_ss)
f_mw_ss2 = (complex double*)calloc((L+1)*(2*L-1), sizeof(complex double));
SSHT_ERROR_MEM_ALLOC_CHECK(f_mw_ss2)
//not freed currenctly


  // Write program name.
  printf("\n");
  printf("SSHT test program (C implementation)\n");
  printf("================================================================\n");

  // Run algorithm error and timing tests.
  for (irepeat = 0; irepeat < NREPEAT; irepeat++) {

    // If spin=0 run tests on algorithms optimised for real spin=0 signal.
    if (spin == 0) {

      // =========================================================================
      // MW real spin=0
      printf("MW real test no. %d\n", irepeat);

      ssht_test_gen_flm_real(flm_orig, L, seed);
      time_start = clock();
      ssht_core_mw_inverse_sov_sym_real(f_mw_real, flm_orig, L, verbosity);      
      time_end = clock();
      durations_inverse_mw_real[irepeat] = 
	(time_end - time_start) / (double)CLOCKS_PER_SEC;
      
      time_start = clock();
      ssht_core_mw_forward_sov_conv_sym_real(flm_syn, f_mw_real, L, verbosity);      
      time_end = clock();
      durations_forward_mw_real[irepeat] = 
	(time_end - time_start) / (double)CLOCKS_PER_SEC;

      errors_mw_real[irepeat] = 0.0;
      for (i = 0; i < L*L; i++)
	errors_mw_real[irepeat] = 
	  MAX(cabs(flm_orig[i] - flm_syn[i]), errors_mw_real[irepeat]);

      printf(" duration_inverse (s) = %40.4f\n", 
	     durations_inverse_mw_real[irepeat]);
      printf(" duration_forward (s) = %40.4f\n", 
	     durations_forward_mw_real[irepeat]);
      printf(" error                = %40.5e\n\n", 
	     errors_mw_real[irepeat]);

      // =========================================================================
      // GL real spin=0
      printf("GL real test no. %d\n", irepeat);

      ssht_test_gen_flm_real(flm_orig, L, seed);
      time_start = clock();
      ssht_core_gl_inverse_sov_real(f_gl_real, flm_orig, L, verbosity);      
      time_end = clock();
      durations_inverse_gl_real[irepeat] = 
	(time_end - time_start) / (double)CLOCKS_PER_SEC;

      time_start = clock();      
      ssht_core_gl_forward_sov_real(flm_syn, f_gl_real, L, verbosity);      
      time_end = clock();
      durations_forward_gl_real[irepeat] = 
	(time_end - time_start) / (double)CLOCKS_PER_SEC;

      errors_gl_real[irepeat] = 0.0;
      for (i = 0; i < L*L; i++)
	errors_gl_real[irepeat] = 
	  MAX(cabs(flm_orig[i] - flm_syn[i]), errors_gl_real[irepeat]);

      printf(" duration_inverse (s) = %40.4f\n", 
	     durations_inverse_gl_real[irepeat]);
      printf(" duration_forward (s) = %40.4f\n", 
	     durations_forward_gl_real[irepeat]);
      printf(" error                = %40.5e\n\n", 
	     errors_gl_real[irepeat]);

      // =========================================================================
      // DH real spin=0
      printf("DH real test no. %d\n", irepeat);

      ssht_test_gen_flm_real(flm_orig, L, seed);
      time_start = clock();
      ssht_core_dh_inverse_sov_real(f_dh_real, flm_orig, L, verbosity);
      time_end = clock();
      durations_inverse_dh_real[irepeat] = 
	(time_end - time_start) / (double)CLOCKS_PER_SEC;
      
      time_start = clock();
      ssht_core_dh_forward_sov_real(flm_syn, f_dh_real, L, verbosity);
      time_end = clock();
      durations_forward_dh_real[irepeat] = 
	(time_end - time_start) / (double)CLOCKS_PER_SEC;

      errors_dh_real[irepeat] = 0.0;
      for (i = 0; i < L*L; i++)
	errors_dh_real[irepeat] = 
	  MAX(cabs(flm_orig[i] - flm_syn[i]), errors_dh_real[irepeat]);

      printf(" duration_inverse (s) = %40.4f\n", 
	     durations_inverse_dh_real[irepeat]);
      printf(" duration_forward (s) = %40.4f\n", 
	     durations_forward_dh_real[irepeat]);
      printf(" error                = %40.5e\n\n", 
	     errors_dh_real[irepeat]);

    }

    // =========================================================================
    // MW
    printf("MW test no. %d\n", irepeat);

    ssht_test_gen_flm_complex(flm_orig, L, spin, seed);
    time_start = clock();
    //ssht_core_mw_inverse_sov_sym(f_mw, flm_orig, L, spin, verbosity);    

    ssht_core_mw_inverse_sov_sym_ss(f_mw_ss, flm_orig, L, spin, verbosity);    
    ssht_core_mwdirect_inverse_ss(f_mw_ss2, flm_orig, L, spin, verbosity);    

    time_end = clock();
    durations_inverse_mw[irepeat] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

    time_start = clock();
    ssht_core_mw_forward_sov_conv_sym(flm_syn, f_mw, L, spin, verbosity);
    time_end = clock();
    durations_forward_mw[irepeat] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

    errors_mw[irepeat] = 0.0;
    for (i = 0; i < L*L; i++)
	errors_mw[irepeat] = 
	  MAX(cabs(flm_orig[i] - flm_syn[i]), errors_mw[irepeat]);

    printf(" duration_inverse (s) = %40.4f\n", 
	   durations_inverse_mw[irepeat]);
    printf(" duration_forward (s) = %40.4f\n", 
	   durations_forward_mw[irepeat]);
    printf(" error                = %40.5e\n\n", 
	   errors_mw[irepeat]);

    // =========================================================================
    // GL
    printf("GL test no. %d\n", irepeat);

    ssht_test_gen_flm_complex(flm_orig, L, spin, seed);
    time_start = clock();
    ssht_core_gl_inverse_sov(f_gl, flm_orig, L, spin, verbosity);
    time_end = clock();
    durations_inverse_gl[irepeat] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

    time_start = clock();
    ssht_core_gl_forward_sov(flm_syn, f_gl, L, spin, verbosity);
    time_end = clock();
    durations_forward_gl[irepeat] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

    errors_gl[irepeat] = 0.0;
    for (i = 0; i < L*L; i++)
	errors_gl[irepeat] = 
	  MAX(cabs(flm_orig[i] - flm_syn[i]), errors_gl[irepeat]);

    printf(" duration_inverse (s) = %40.4f\n", 
	   durations_inverse_gl[irepeat]);
    printf(" duration_forward (s) = %40.4f\n", 
	   durations_forward_gl[irepeat]);
    printf(" error                = %40.5e\n\n", 
	   errors_gl[irepeat]);

    // =========================================================================
    // DH
    printf("DH test no. %d\n", irepeat);

    ssht_test_gen_flm_complex(flm_orig, L, spin, seed);
    time_start = clock();
    ssht_core_dh_inverse_sov(f_dh, flm_orig, L, spin, verbosity);
    time_end = clock();
    durations_inverse_dh[irepeat] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

    time_start = clock();
    ssht_core_dh_forward_sov(flm_syn, f_dh, L, spin, verbosity);
    time_end = clock();
    durations_forward_dh[irepeat] = (time_end - time_start) / (double)CLOCKS_PER_SEC;

    errors_dh[irepeat] = 0.0;
    for (i = 0; i < L*L; i++)
	errors_dh[irepeat] = 
	  MAX(cabs(flm_orig[i] - flm_syn[i]), errors_dh[irepeat]);

    printf(" duration_inverse (s) = %40.4f\n", 
	   durations_inverse_dh[irepeat]);
    printf(" duration_forward (s) = %40.4f\n", 
	   durations_forward_dh[irepeat]);
    printf(" error                = %40.5e\n\n", 
	   errors_dh[irepeat]);
   
  }

  // =========================================================================
  // Summarise results

  printf("================================================================\n");
  printf("Summary\n\n");
  printf("NREPEAT               = %40d\n", NREPEAT);
  printf("L                     = %40d\n", L);
  printf("spin                  = %40d\n\n", spin);

  if (spin == 0) {

    printf("MW real\n");
    tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_forward_mw_real[i];
    printf(" Average forward transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
    tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_inverse_mw_real[i];
    printf(" Average inverse transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
    tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += errors_mw_real[i];
    printf(" Average max error                  = %26.5e\n\n", tmp/(double)NREPEAT);

    printf("GL real\n");
    tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_forward_gl_real[i];
    printf(" Average forward transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
    tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_inverse_gl_real[i];
    printf(" Average inverse transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
    tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += errors_gl_real[i];
    printf(" Average max error                  = %26.5e\n\n", tmp/(double)NREPEAT);

    printf("DH real\n");
    tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_forward_dh_real[i];
    printf(" Average forward transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
    tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_inverse_dh_real[i];
    printf(" Average inverse transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
    tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += errors_dh_real[i];
    printf(" Average max error                  = %26.5e\n\n", tmp/(double)NREPEAT);

  }

  printf("MW\n");
  tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_forward_mw[i];
  printf(" Average forward transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
  tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_inverse_mw[i];
  printf(" Average inverse transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
  tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += errors_mw[i];
  printf(" Average max error                  = %26.5e\n\n", tmp/(double)NREPEAT);

  printf("GL\n");
  tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_forward_gl[i];
  printf(" Average forward transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
  tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_inverse_gl[i];
  printf(" Average inverse transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
  tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += errors_gl[i];
  printf(" Average max error                  = %26.5e\n\n", tmp/(double)NREPEAT);

  printf("DH\n");
  tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_forward_dh[i];
  printf(" Average forward transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
  tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += durations_inverse_dh[i];
  printf(" Average inverse transform time (s) = %26.4f\n", tmp/(double)NREPEAT);
  tmp = 0.0; for (i=0; i<NREPEAT; i++) tmp += errors_dh[i];
  printf(" Average max error                  = %26.5e\n\n", tmp/(double)NREPEAT);

  // Free memory.
  free(flm_orig);
  free(flm_syn);
  free(f_mw);
  free(f_mw_real);
  free(f_gl);
  free(f_gl_real);
  free(f_dh);
  free(f_dh_real);

  return 0;
}


/*!  
 * Generate random spherical harmonic coefficients of a real spin=0
 * signal.
 *
 * \param[out] flm Random spherical harmonic coefficients generated.
 * \param[in] L Harmonic band-limit.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_test_gen_flm_real(complex double *flm, int L, int seed) {

  int el, m, msign, i, i_op;

  for (el=0; el<L; el++) {
    m = 0;
    ssht_sampling_elm2ind(&i, el, m);
    flm[i] = (2.0*ran2_dp(seed) - 1.0);
    for (m=1; m<=el; m++) {
      ssht_sampling_elm2ind(&i, el, m);
      flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);
      ssht_sampling_elm2ind(&i_op, el, -m);
      msign = m & 1;
      msign = 1 - msign - msign; // (-1)^m
      flm[i_op] = msign * conj(flm[i]);
    }
  }

}


/*!  
 * Generate random spherical harmonic coefficients of a complex
 * signal.
 *
 * \param[out] flm Random spherical harmonic coefficients generated.
 * \param[in] L Harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_test_gen_flm_complex(complex double *flm, int L, int spin, int seed) {

  int i, i_lo;

  ssht_sampling_elm2ind(&i_lo, abs(spin), 0);
  for (i=i_lo; i<L*L; i++) 
    flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);

}


/*!  
 * Generate uniform deviate in range [0,1) given seed. (Using double
 * precision.)
 *
 * \notes Uniform deviate (Num rec 1992, chap 7.1), original routine
 * said to be 'perfect'.
 *
 * \param[in] idum Seed.
 * \retval ran_dp Generated uniform deviate.
 *
 * \author Jason McEwen
 */
double ran2_dp(int idum) {

  int IM1=2147483563,IM2=2147483399,IMM1=IM1-1, 
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, 
    NTAB=32,NDIV=1+IMM1/NTAB;

  double AM=1./IM1,EPS=1.2e-7,RNMX=1.-EPS;
  int j,k;
  static int iv[32],iy,idum2 = 123456789; 
  // N.B. in C static variables are initialised to 0 by default.

  if (idum <= 0) {
    idum= (-idum>1 ? -idum : 1); // max(-idum,1);
    idum2=idum;
    for(j=NTAB+8;j>=1;j--) {
      k=idum/IQ1;
      idum=IA1*(idum-k*IQ1)-k*IR1;
      if (idum < 0) idum=idum+IM1;
      if (j < NTAB) iv[j-1]=idum;
    }
    iy=iv[0];
  }
  k=idum/IQ1;
  idum=IA1*(idum-k*IQ1)-k*IR1;
  if (idum < 0) idum=idum+IM1;
  k=idum2/IQ2;
  idum2=IA2*(idum2-k*IQ2)-k*IR2;
  if (idum2 < 0) idum2=idum2+IM2;
  j=1+iy/NDIV;
  iy=iv[j-1]-idum2;
  iv[j-1]=idum;
  if(iy < 1)iy=iy+IMM1;
  return (AM*iy < RNMX ? AM*iy : RNMX); // min(AM*iy,RNMX);

}
