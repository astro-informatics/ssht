


#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <complex.h>
#include <fftw3.h>
// Since complex.h included before fftw3.h, fftw_complex is
// defined to be the native complex type in complex.h.
#include "ssht_sampling.h"
#include "ssht_error.h"


double ran2_dp(int idum);


int main(int argc, char *argv[]) {

  complex *flm;
  int L = 128;
  int i;
  int seed = 100;

  printf("hello world!\n");





  flm = (complex*)calloc(L*L, sizeof(complex));
  SSHT_ERROR_MEM_ALLOC_CHECK(flm)

 
  for(i=0; i<10; i++) {
    printf("ran2_dp = %f\n", ran2_dp(seed));
  }

  // Free memory.
  free(flm);
}


void ssht_test_gen_flm_complex(complex *flm, int L, int spin, int seed) {

  int i, i_lo;

  ssht_sampling_elm2ind(&i_lo, abs(spin), 0);

  for (i=i_lo; i<L*L; i++) {
    flm[i] = 1.0 + I * 1.0;
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
