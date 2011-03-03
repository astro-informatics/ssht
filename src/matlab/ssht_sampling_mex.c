// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#include <ssht.h>
#include "ssht_mex.h"
#include <string.h>
#include "mex.h"


/**
 * Compute sampling positions.
 *
 * Usage: 
 *   [n, ntheta, nphi, thetas, phis] = ssht_sampling_mex(L, method)
 *
 * \author Jason McEwen
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int L, ntheta, nphi, n, t, p, iin;
  double *thetas, *phis, *weights_unused;
  int len, iout = 0;
  char method[SSHT_STRING_LEN];

  /* Check number of arguments. */
  if(nrhs!=2) {
    mexErrMsgIdAndTxt("ssht_sampling_mex:InvalidInput:nrhs",
		      "Require two inputs.");
  }
  if(nlhs!=5) {
    mexErrMsgIdAndTxt("ssht_sampling_mex:InvalidOutput:nlhs",
		      "Require five outputs.");
  }

  /* Parse harmonic band-limit L. */
  iin = 0;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_sampling_mex:InvalidInput:bandLimit",
		      "Harmonic band-limit must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0) {
    mexErrMsgIdAndTxt("ssht_sampling_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit must be positive integer.");
  }

  /* Parse method. */
  iin = 1;
  if( !mxIsChar(prhs[iin]) ) {
    mexErrMsgIdAndTxt("ssht_sampling_mex:InvalidInput:methodChar",
		      "Method must be string.");
  }
  len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
  if (len >= SSHT_STRING_LEN) 
    mexErrMsgIdAndTxt("ssht_sampling_mex:InvalidInput:methodTooLong",
		      "Method exceeds string length.");
  mxGetString(prhs[iin], method, len);

  /* Compute sample positions. */
  if (strcmp(method, SSHT_SAMPLING_MW) == 0) {

    ntheta = ssht_sampling_mw_ntheta(L);
    nphi = ssht_sampling_mw_nphi(L);
    n = ssht_sampling_mw_n(L);
    plhs[iout++] = mxCreateDoubleScalar(n);
    plhs[iout++] = mxCreateDoubleScalar(ntheta);
    plhs[iout++] = mxCreateDoubleScalar(nphi);
      
    plhs[iout] = mxCreateDoubleMatrix(ntheta, 1, mxREAL);
    thetas = mxGetPr(plhs[iout++]);
    for (t=0; t<ntheta; t++)
      thetas[t] = ssht_sampling_mw_t2theta(t, L);
    plhs[iout] = mxCreateDoubleMatrix(nphi, 1, mxREAL);
    phis = mxGetPr(plhs[iout++]);
    for (p=0; p<nphi; p++)
      phis[p] = ssht_sampling_mw_p2phi(p, L);

  }
  else if (strcmp(method, SSHT_SAMPLING_MWSS) == 0) {
      
    ntheta = ssht_sampling_mw_ss_ntheta(L);
    nphi = ssht_sampling_mw_ss_nphi(L);
    n = ssht_sampling_mw_ss_n(L);
    plhs[iout++] = mxCreateDoubleScalar(n);
    plhs[iout++] = mxCreateDoubleScalar(ntheta);
    plhs[iout++] = mxCreateDoubleScalar(nphi);

    plhs[iout] = mxCreateDoubleMatrix(ntheta, 1, mxREAL);
    thetas = mxGetPr(plhs[iout++]);
    for (t=0; t<ntheta; t++)
      thetas[t] = ssht_sampling_mw_ss_t2theta(t, L);
    plhs[iout] = mxCreateDoubleMatrix(nphi, 1, mxREAL);
    phis = mxGetPr(plhs[iout++]);
    for (p=0; p<nphi; p++)
      phis[p] = ssht_sampling_mw_ss_p2phi(p, L);

  }
  else if (strcmp(method, SSHT_SAMPLING_GL) == 0) {

    ntheta = ssht_sampling_gl_ntheta(L);
    nphi = ssht_sampling_gl_nphi(L);
    n = ssht_sampling_gl_n(L);
    plhs[iout++] = mxCreateDoubleScalar(n);
    plhs[iout++] = mxCreateDoubleScalar(ntheta);
    plhs[iout++] = mxCreateDoubleScalar(nphi);

    plhs[iout] = mxCreateDoubleMatrix(ntheta, 1, mxREAL);
    thetas = mxGetPr(plhs[iout++]);
    weights_unused = (double*)malloc(L*sizeof(double));
    ssht_sampling_gl_thetas_weights(thetas, weights_unused, L);
    free(weights_unused);
    plhs[iout] = mxCreateDoubleMatrix(nphi, 1, mxREAL);
    phis = mxGetPr(plhs[iout++]);
    for (p=0; p<nphi; p++)
      phis[p] = ssht_sampling_gl_p2phi(p, L);

  }
  else if (strcmp(method, SSHT_SAMPLING_DH) == 0) {

    ntheta = ssht_sampling_dh_ntheta(L);
    nphi = ssht_sampling_dh_nphi(L);
    n = ssht_sampling_dh_n(L);
    plhs[iout++] = mxCreateDoubleScalar(n);
    plhs[iout++] = mxCreateDoubleScalar(ntheta);
    plhs[iout++] = mxCreateDoubleScalar(nphi);

    plhs[iout] = mxCreateDoubleMatrix(ntheta, 1, mxREAL);
    thetas = mxGetPr(plhs[iout++]);
    for (t=0; t<ntheta; t++)
      thetas[t] = ssht_sampling_dh_t2theta(t, L);
    plhs[iout] = mxCreateDoubleMatrix(nphi, 1, mxREAL);
    phis = mxGetPr(plhs[iout++]);
    for (p=0; p<nphi; p++)
      phis[p] = ssht_sampling_dh_p2phi(p, L);

  }
  else {

    mexErrMsgIdAndTxt("ssht_sampling_mex:InvalidInput:method",
		      "Method invalid.");

  } 

}
