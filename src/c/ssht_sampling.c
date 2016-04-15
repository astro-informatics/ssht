// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details

/*! 
 * \file ssht_sampling.c
 * Functionality to define sample positions for various algorithms,
 * to compute weights and to convert 1D and 2D harmonic indices.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */


#include <complex.h>
#include <math.h>
#include "ssht_types.h"

#define ABS(x) ((x) >= 0 ? (x) : -(x))

void gauleg(double x1, double x2, double *x, double *w, int n);


//============================================================================
// Sampling weights
//============================================================================


/*!
 * Compute weights for toroidal extension.
 *
 * \param[in]  p Integer index to compute weight for.
 * \param[out] w Corresponding weight.
 * \retval Weight value computed.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
complex double ssht_sampling_weight_mw(int p) {

  if (p == 1) {
    return I * SSHT_PION2;
  }
  else if (p == -1) {
     return - I * SSHT_PION2;
  }
  else if (p % 2 == 0) {
    return 2.0 / (1.0 - p*p);
  }
  else {
    return 0.0;
  }
  
}
  

/*!
 * Compute Driscoll and Healy weights.
 *
 * \param[in] theta_t Theta value to compute weight for.
 * \param[in] L Harmonic band-limit.
 * return w Weight value computed.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

double ssht_sampling_weight_dh(double theta_t, int L) {

  double w;
  int k;

  w = 0.0;
  for (k=0; k<=L-1; k++)
    w += sin((2.0*k+1.0)*theta_t) / (double)(2.0*k+1.0);

  w *= 2.0/((double) L) * sin(theta_t);

  return w;

}


/*!
 * Compute Gauss-Legendre theta positions (roots of Legendre
 * polynomials) and corresponding weights.
 *
 * \param[out] thetas L theta positions (memory must already be
 * allocated to store the L theta positions).
 * \param[out] weights Corresponding weights (memory must already be
 * allocated to store the L weights corresponding to each theta position).
 * \param[in] L Harmonic band-limit.
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_sampling_gl_thetas_weights(double *thetas, double *weights, int L) {
   
  int t;
  double tmp;

  gauleg(-1.0, 1.0, thetas, weights, L);

  for (t=0; t<=L-1; t++)
    thetas[t] = acos(thetas[t]);

  // Reorder thetas array.
  for (t=0; t<=(L-1)/2; t++) {
    tmp = thetas[t];
    thetas[t] = thetas[L-1-t];
    thetas[L-1-t] = tmp;
  }

}


/*!
 * Given the lower and upper limits of integration x1 and x2, this
 * routine returns arrays x[1..n] and w[1..n] of length n,
 * containing the abscissas and weights of the Gauss-Legendre
 * n-point quadrature formula.
 *
 * \param[in] x1 Lower bound of range.
 * \param[in] x2 Upper bound of range.
 * \param[out] x Node positions (i.e. roots of Legendre polynomials).
 * \param[out] w Corresponding weights.
 * \param[in] n Number of points (note memory must already be
 * allocated for x and w to store n terms).
 *
 * \author Numerical recipes.
 */
void gauleg(double x1, double x2, double *x, double *w, int n) {

  double EPS = 1e-14;
  int i,j,m;
  double p1,p2,p3,pp,xl,xm,z,z1;

  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);

  for (i=1; i<=m; i++) {
    z=cos(3.141592654*(i-.25)/(n+.5));
    do {
      p1=1.0;
      p2=0.0;
      for (j=1; j<=n; j++) {
	p3=p2;
	p2=p1;
	p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
      }
      pp=n*(z*p1-p2)/(z*z-1.0);
      z1=z;
      z=z1-p1/pp;
    } while (ABS(z-z1) > EPS);
    x[i-1]=xm-xl*z;
    x[n+1-i-1]=xm+xl*z;
    w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
    w[n+1-i-1]=w[i-1];
  }

}


//============================================================================
// Sampling relations
//============================================================================


/*!
 * Compute number of theta samples for McEwen and Wiaux sampling.
 * 
 * /note Computes number of samples in (0,pi], *not* over extended
 * domain.
 *
 * \param[in] L Harmonic band-limit.
 * \retval ntheta Number of theta samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_mw_ntheta(int L) {

  return L;

}


/*!
 * Convert theta index to angle for McEwen and Wiaux sampling.
 *
 * \note
 *  - t ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
 *
 * \param[in] t Theta index.
 * \param[in] L Harmonic band-limit.
 * \retval theta Theta angle.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double ssht_sampling_mw_t2theta(int t, int L) {

  return (2.0*t + 1.0) * SSHT_PI / (2.0*L - 1.0);

}


/*!
 * Compute number of phi samples for McEwen and Wiaux sampling.
 *
 * \param[in] L Harmonic band-limit.
 * \retval nphi Number of phi samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_mw_nphi(int L) {

  return 2*L-1;

}


/*!
 * Convert phi index to angle for McEwen and Wiaux sampling.
 *
 * \note
 *  - p ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
 *
 * \param[in] p Phi index.
 * \param[in] L Harmonic band-limit.
 * \retval phi Phi angle.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double ssht_sampling_mw_p2phi(int p, int L) {

  return 2.0 * p * SSHT_PI / (2.0*L - 1.0);

}


/*!
 * Compute total number of samples for McEwen and Wiaux sampling.
 * 
 * /note Computes number of samples on sphere, *not* over extended
 * domain.
 *
 * \param[in] L Harmonic band-limit.
 * \retval n Number of samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_mw_n(int L) {

  return (L-1)*(2*L-1) + 1;

}


/*!
 * Convert theta index to angle for McEwen and Wiaux symmetric
 * sampling.
 *
 * \note
 *  - t ranges from [0 .. 2*L-1] => 2*L points in (0,2*pi).
 *
 * \param[in] t Theta index.
 * \param[in] L Harmonic band-limit.
 * \retval theta Theta angle.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double ssht_sampling_mw_ss_t2theta(int t, int L) {

  return 2.0 * t * SSHT_PI / (2.0 * L);

}


/*!
 * Compute number of theta samples for McEwen and Wiaux symmetric
 * sampling.
 * 
 * /note Computes number of samples in [0,pi], *not* over extended
 * domain.
 *
 * \param[in] L Harmonic band-limit.
 * \retval ntheta Number of theta samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_mw_ss_ntheta(int L) {

  return L+1;

}


/*!
 * Convert phi index to angle for McEwen and Wiaux symmetric sampling.
 *
 * \note
 *  - p ranges from [0 .. 2*L-1] => 2*L points in [0,2*pi).
 *
 * \param[in] p Phi index.
 * \param[in] L Harmonic band-limit.
 * \retval phi Phi angle.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double ssht_sampling_mw_ss_p2phi(int p, int L) {

  return 2.0 * p * SSHT_PI / (2.0*L);

}


/*!
 * Compute number of phi samples for McEwen and Wiaux symmetric
 * sampling.
 *
 * \param[in] L Harmonic band-limit.
 * \retval nphi Number of phi samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_mw_ss_nphi(int L) {

  return 2*L;

}


/*!
 * Compute total number of samples for McEwen and Wiaux symmetric
 * sampling.
 * 
 * /note Computes number of samples on sphere, *not* over extended
 * domain.
 *
 * \param[in] L Harmonic band-limit.
 * \retval n Number of samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_mw_ss_n(int L) {

  return (L-1)*(2*L) + 2;

}


/*!
 * Convert theta index to angle for Driscoll and Healy sampling.
 *
 * \note
 *  - t ranges from [0 .. 2*L-1] => 2*L points in (0,pi).
 *
 * \param[in] t Theta index.
 * \param[in] L Harmonic band-limit.
 * \retval theta Theta angle.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double ssht_sampling_dh_t2theta(int t, int L) {

  return (2.0*t + 1.0) * SSHT_PI / (4.0*L);

}


/*!
 * Compute number of theta samples for Driscoll and Healy sampling.
 *
 * \param[in] L Harmonic band-limit.
 * \retval ntheta Number of theta samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_dh_ntheta(int L) {

  return 2*L;

}


/*!
 * Convert phi index to angle for Driscoll and Healy sampling.
 *
 * \note
 *  - p ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
 *
 * \param[in] p Phi index.
 * \param[in] L Harmonic band-limit.
 * \retval phi Phi angle.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double ssht_sampling_dh_p2phi(int p, int L) {

  return 2.0 * p * SSHT_PI / (2.0*L - 1.0);

}


/*!
 * Compute number of phi samples for Driscoll and Healy sampling.
 *
 * \param[in] L Harmonic band-limit.
 * \retval nphi Number of phi samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_dh_nphi(int L) {

  return 2*L-1;

}


/*!
 * Compute total number of samples for Driscoll and Healy sampling.
 * 
 * \param[in] L Harmonic band-limit.
 * \retval n Number of samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_dh_n(int L) {

  return 2*L*(2*L-1);

}


/*!
 * Compute number of theta samples for Gauss-Legendre sampling.
 *
 * \param[in] L Harmonic band-limit.
 * \retval ntheta Number of theta samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_gl_ntheta(int L) {

  return L;

}


/*!
 * Convert phi index to angle for Gauss-Legendre sampling.
 *
 * \note
 *  - p ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
 *
 * \param[in] p Phi index.
 * \param[in] L Harmonic band-limit.
 * \retval phi Phi angle.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double ssht_sampling_gl_p2phi(int p, int L) {

  return 2.0 * p * SSHT_PI / (2.0*L - 1.0);

}


/*!
 * Compute number of phi samples for Gauss-Legendre sampling.
 *
 * \param[in] L Harmonic band-limit.
 * \retval nphi Number of phi samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_gl_nphi(int L) {

  return 2*L-1;

}


/*!
 * Compute total number of samples for Gauss-Legendre sampling.
 *
 * \param[in] L Harmonic band-limit.
 * \retval n Number of samples.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_sampling_gl_n(int L) {

  return L*(2*L-1);

}


//============================================================================
// Harmonic index relations
//============================================================================


/*!
 * Convert (el,m) harmonic indices to 1D index used to access flm
 * array.
 *
 * \note Index ranges are as follows:  
 *  - el ranges from [0 .. L-1].
 *  - m ranges from [-el .. el].
 *  - ind ranges from [0 .. L**2-1].
 *
 * \param[out] ind 1D index to access flm array [output].
 * \param[in]  el  Harmonic index [input].
 * \param[in]  m   Azimuthal harmonic index [input].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_sampling_elm2ind(int *ind, int el, int m) {

  *ind = el * el + el + m;

}


/*!
 * Convert 1D index used to access flm array to (el,m) harmonic
 * indices.
 *
 * \note Index ranges are as follows:  
 *  - el ranges from [0 .. L-1].
 *  - m ranges from [-el .. el].
 *  - ind ranges from [0 .. L**2-1].
 *
 * \param[in]  ind 1D index to access flm array [output].
 * \param[out] el  Harmonic index [input].
 * \param[out] m   Azimuthal harmonic index [input].
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_sampling_ind2elm(int *el, int *m, int ind) {

  *el = sqrt(ind);
  *m = ind - (*el)*(*el) - (*el);

}

     
