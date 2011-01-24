


#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <complex.h>
#include <fftw3.h>
// Since complex.h included before fftw3.h, fftw_complex is
// defined to be the native complex type in complex.h.
#include "ssht_sampling.h"
#include "ssht_error.h"

#define NREPEAT 2

double ran2_dp(int idum);
void ssht_test_gen_flm_complex(complex double *flm, int L, int spin, int seed);
void ssht_test_gen_flm_real(complex double *flm, int L, int seed);


int main(int argc, char *argv[]) {

  complex double *flm_orig, *flm_syn;
  complex double *f_mw;
  int L = 128;
  int spin = 0;
  int irepeat;
  int seed = 1;
  int verbosity = 0;

  double max_err[NREPEAT];
  double tmp;




  int i;

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

  // Write program name.
  printf("\n");
  printf("SSHT test program (C implementation)\n");
  printf("===============================================================\n");

  // Run algorithm error and timing tests.
  for (irepeat = 0; irepeat < NREPEAT; irepeat++) {

    // If spin=0 run tests on algorithms optimised for real spin=0 signal.
    if (spin == 0) {

      printf("TODO\n");
      printf("L=%d, spin=%d\n", L, spin);

    }

    // =========================================================================
    // MW
    printf("MW test no. %d\n", irepeat);

    ssht_test_gen_flm_complex(flm_orig, L, spin, seed);

    ssht_core_mw_inverse_sov_sym(f_mw, flm_orig, L, spin, verbosity);
    ssht_core_mw_forward_sov_conv_sym(flm_syn, f_mw, L, spin, verbosity);

    max_err[irepeat] = 0.0;
    for (i = 0; i < L*L; i++) {
      tmp = cabs(flm_orig[i] - flm_syn[i]);
      max_err[irepeat] = tmp > max_err[irepeat] ? tmp : max_err[irepeat];
    }
    printf(" error %40.5e\n", max_err[irepeat]);


  }







  for(i=0; i<10; i++) {
    printf("ran2_dp = %f\n", ran2_dp(seed));
  }

  // Free memory.
  free(flm_orig);
  free(flm_syn);
  free(f_mw);
}






void ssht_test_gen_flm_complex(complex double *flm, int L, int spin, int seed) {

  int i, i_lo;

  ssht_sampling_elm2ind(&i_lo, abs(spin), 0);
  for (i=i_lo; i<L*L; i++) 
    flm[i] = (2.0*ran2_dp(seed) - 1.0) + I * (2.0*ran2_dp(seed) - 1.0);

}



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



// update to use ran c generator from NR

double ran2_dp(int idum) {

  //int IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV;
  //double AM,EPS,RNMX;

// this is c++ not c, i.e. use of const!
  int IM1=2147483563,IM2=2147483399,IMM1=IM1-1, 
    IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, 
    NTAB=32,NDIV=1+IMM1/NTAB;

  const double AM=1./IM1,EPS=1.2e-7,RNMX=1.-EPS;

  int j,k;
  static int iv[32],iy,idum2 = 123456789; 
// iv[NTAB]
// static variables initialised to 0 by default

  if (idum <= 0) {
    idum= (-idum>1 ? -idum : 1);//max(-idum,1);
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
  return (AM*iy < RNMX ? AM*iy : RNMX); //min(AM*iy,RNMX);

}
