// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2013  Jason McEwen
// See LICENSE.txt for license details


#include <ssht.h>
#include "ssht_mex.h"
#include <string.h>
#include "mex.h"


/**
 * Compute Wigner small-d functions using Risbo's method.
 *
 * Usage: 
 *   dl = ssht_dl_mex(dl, L, el, theta, sqrt_tbl, signs);
 *
 * \author Jason McEwen
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int iin = 0, iout = 0, i;
  double *dl = NULL, *dl_out = NULL;
  int dl_m, dl_n;
  int L, el;
  double theta;
  double *sqrt_tbl, *signs;
  int sqrt_tbl_m, sqrt_tbl_n, signs_m, signs_n;

  /* Check number of arguments. */
  if(nrhs!=6) {
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:nrhs",
		      "Require 6 inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidOutput:nlhs",
		      "Require 1 output.");
  }

  /* Parse dl array. */
  iin = 0;
  if( !mxIsDouble(prhs[iin]) ) {
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:dl",
		      "Wigner function values must be doubles.");
  }
  dl_m = mxGetM(prhs[iin]);
  dl_n = mxGetN(prhs[iin]);  
  dl = mxGetPr(prhs[iin]); 

  /* Parse harmonic band-limit L. */
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:bandLimit",
		      "Harmonic band-limit must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit must be positive integer.");
  if (dl_m != 2*L-1 || dl_n != 2*L-1)
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:dl",
		      "Invalid size for dl matrix.");

  /* Parse harmonic degree el. */
  iin = 2;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:el",
		      "Harmonic degree el must be integer.");
  }
  el = (int)mxGetScalar(prhs[iin]);
  if (el >= L || el < 0)
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:elInvalid",
		      "Harmonic degree el invalid.");

  /* Parse theta. */
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:theta",
		      "Theta must be real double.");
  }
  theta = mxGetScalar(prhs[iin]);
  if (theta < 0 || theta > SSHT_PI)
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:thetaInvalid",
		      "Harmonic band-limit must be positive integer.");

  /* Parse sqrt_tbl array. */
  iin = 4;
  sqrt_tbl_m = mxGetM(prhs[iin]);
  sqrt_tbl_n = mxGetN(prhs[iin]);  
  if( !mxIsDouble(prhs[iin]) ||
      sqrt_tbl_m * sqrt_tbl_n != 2*(L-1)+2
      || (sqrt_tbl_m != 1 && sqrt_tbl_n != 1) ) {
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:sqrt_tbl",
		      "Square-root table must be double array.");
  }
  sqrt_tbl = mxGetPr(prhs[iin]); 

  /* Parse signs array. */
  iin = 5;
  signs_m = mxGetM(prhs[iin]);
  signs_n = mxGetN(prhs[iin]);  
  if( signs_m * signs_n != L+1 
      || (signs_m != 1 && signs_n != 1)) {
    mexErrMsgIdAndTxt("ssht_dl_mex:InvalidInput:signs",
		      "Sign table must be array.");
  }
  signs = mxGetPr(prhs[iin]); 

  /* Compute Wigner small-d functions. */
  ssht_dl_beta_risbo_half_table(dl, theta, L, SSHT_DL_FULL, 
		el, sqrt_tbl, signs);

  /* Copy dl to output. */
  plhs[iout] = mxCreateDoubleMatrix(2*L-1, 2*L-1, mxREAL);
  dl_out = mxGetPr(plhs[iout]);
  for (i=0; i<(2*L-1)*(2*L-1); i++)
    dl_out[i] = dl[i];

}
