#include "ssht.h"
#include "utilities.h"

#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

/*!
 * Test for null vector
 *
 * \param[in]  X vector of complex double
 * \param[in]  n length of X.
 * \retval     Y returns int 0 if all zeros, 1 if non-zero
 *             i.e. 0 = fail, 1 = pass.
 */
int null_test(const complex double *X, int n) {
  int Y = 0;
  for (int i = 0; i < n; ++i) {
    if (cabs(X[i]) != 0.0) {
      Y = 1;
      i = n;
    }
  }
  return Y;
}

/*!
 * Test for nan vector
 *
 * \param[in]  X vector of complex double
 * \param[in]  n length of X.
 * \retval     Y returns int 1 if no nans, 0 if nan entires exist.
 *             i.e. 0 = fail, 1 = pass.
 */
int nan_test(const complex double *X, int n) {
  int Y = 1;
  for (int i = 0; i < n; ++i) {
    if (X[i] != X[i]) {
      Y = 0;
      i = n;
    }
  }
  return Y;
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
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void gen_flm_real(complex double *flm, int L, int seed) {

  int el, m, msign, i, i_op;

  for (el = 0; el < L; el++) {
    m = 0;
    ssht_sampling_elm2ind(&i, el, m);
    flm[i] = (2.0 * ran2_dp(seed) - 1.0);
    for (m = 1; m <= el; m++) {
      ssht_sampling_elm2ind(&i, el, m);
      flm[i] = (2.0 * ran2_dp(seed) - 1.0) + I * (2.0 * ran2_dp(seed) - 1.0);
      ssht_sampling_elm2ind(&i_op, el, -m);
      msign = m & 1;
      msign = 1 - msign - msign; // (-1)^m
      flm[i_op] = msign * conj(flm[i]);
    }
  }
}

/*!
 * Generate random spherical harmonic coefficients of a real spin=0
 * signal with lower band-limit.
 *
 * \param[out] flm Random spherical harmonic coefficients generated.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void gen_lb_flm_real(complex double *flm, int L0, int L, int seed) {

  int el, m, msign, i, i_op;

  for (el = 0; el < L0; el++) {
    m = 0;
    ssht_sampling_elm2ind(&i, el, m);
    flm[i] = 0.0;
    for (m = 1; m <= el; m++) {
      ssht_sampling_elm2ind(&i, el, m);
      flm[i] = 0.0;
      ssht_sampling_elm2ind(&i, el, -m);
      flm[i] = 0.0;
    }
  }

  for (el = L0; el < L; el++) {
    m = 0;
    ssht_sampling_elm2ind(&i, el, m);
    flm[i] = (2.0 * ran2_dp(seed) - 1.0);
    for (m = 1; m <= el; m++) {
      ssht_sampling_elm2ind(&i, el, m);
      flm[i] = (2.0 * ran2_dp(seed) - 1.0) + I * (2.0 * ran2_dp(seed) - 1.0);
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
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void gen_flm_complex(complex double *flm, int L, int spin, int seed) {

  int i, i_lo;

  ssht_sampling_elm2ind(&i_lo, abs(spin), 0);
  for (i = i_lo; i < L * L; i++)
    flm[i] = (2.0 * ran2_dp(seed) - 1.0) + I * (2.0 * ran2_dp(seed) - 1.0);
}

int max(int a, int b) {
  if (a > b)
    return a;
  return b;
}

/*!
 * Generate random spherical harmonic coefficients of a complex
 * signal with lower band-limit.
 *
 * \param[out] flm Random spherical harmonic coefficients generated.
 * \param[in] L0 Lower harmonic band-limit.
 * \param[in] L Upper harmonic band-limit.
 * \param[in] spin Spin number.
 * \param[in] seed Integer seed required for random number generator.
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void gen_lb_flm_complex(complex double *flm, int L0, int L, int spin, int seed) {

  int i, i_lo;

  ssht_sampling_elm2ind(&i_lo, abs(spin), 0);
  for (i = 0; i < max(i_lo, L0 * L0); i++) {
    flm[i] = 0.0;
  }
  for (i = max(i_lo, L0 * L0); i < L * L; i++) {
    flm[i] = (2.0 * ran2_dp(seed) - 1.0) + I * (2.0 * ran2_dp(seed) - 1.0);
  }
}

/*!
 * Generate uniform deviate in range [0,1) given seed. (Using double
 * precision.)
 *
 * \note Uniform deviate (Num rec 1992, chap 7.1), original routine
 * said to be 'perfect'.
 *
 * \param[in] idum Seed.
 * \retval ran_dp Generated uniform deviate.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double ran2_dp(int idum) {

  int IM1 = 2147483563, IM2 = 2147483399, IMM1 = IM1 - 1, IA1 = 40014, IA2 = 40692,
      IQ1 = 53668, IQ2 = 52774, IR1 = 12211, IR2 = 3791, NTAB = 32,
      NDIV = 1 + IMM1 / NTAB;

  double AM = 1. / IM1, EPS = 1.2e-7, RNMX = 1. - EPS;
  int j, k;
  static int iv[32], iy, idum2 = 123456789;
  // N.B. in C static variables are initialised to 0 by default.

  if (idum <= 0) {
    idum = (-idum > 1 ? -idum : 1); // max(-idum,1);
    idum2 = idum;
    for (j = NTAB + 8; j >= 1; j--) {
      k = idum / IQ1;
      idum = IA1 * (idum - k * IQ1) - k * IR1;
      if (idum < 0)
        idum = idum + IM1;
      if (j < NTAB)
        iv[j - 1] = idum;
    }
    iy = iv[0];
  }
  k = idum / IQ1;
  idum = IA1 * (idum - k * IQ1) - k * IR1;
  if (idum < 0)
    idum = idum + IM1;
  k = idum2 / IQ2;
  idum2 = IA2 * (idum2 - k * IQ2) - k * IR2;
  if (idum2 < 0)
    idum2 = idum2 + IM2;
  j = 1 + iy / NDIV;
  iy = iv[j - 1] - idum2;
  iv[j - 1] = idum;
  if (iy < 1)
    iy = iy + IMM1;
  return (AM * iy < RNMX ? AM * iy : RNMX); // min(AM*iy,RNMX);
}
