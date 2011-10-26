// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#include <ssht.h>
#include "ssht_mex.h"
#include <string.h>
#include "mex.h"


/**
 * Compute forward adjoint transform.
 *
 * Usage: 
 *   [f, f_sp, phi_sp, f_np, phi_np] = ...
 *     ssht_forward_adjoint_mex(flm, L, method, spin, reality, ...
 *                              southPoleExists, northPoleExists);
 *
 * \author Jason McEwen
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  ssht_dl_method_t dl_method = SSHT_DL_TRAPANI;
  int i, L, spin, reality, verbosity=0, flm_m, flm_n;
  int south_pole_exists, north_pole_exists;
  double *flm_real, *flm_imag, *f_real, *f_imag;
  double *fr = NULL;
  complex double *flm = NULL, *f = NULL;
  complex double f_sp = 0.0, f_np = 0.0;
  double fr_sp = 0.0, fr_np = 0.0;
  double phi_sp = 0.0, phi_np = 0.0;
  int ntheta = -1, nphi = -1, t, p;
  int len, iin = 0, iout = 0;
  char method[SSHT_STRING_LEN];

  /* Check number of arguments. */
  if(nrhs!=7) {
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:nrhs",
		      "Require five inputs.");
  }
  if(nlhs!=5) {
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidOutput:nlhs",
		      "Require one output.");
  }

  /* Parse reality. */
  iin = 4;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:reality",
		      "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);

  /* Parse harmonic coefficients flm. */
  iin = 0;
  flm_m = mxGetM(prhs[iin]);
  flm_n = mxGetN(prhs[iin]);
  if (flm_m != 1 && flm_n != 1)
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:flmVector",
		      "Harmonic coefficients must be contained in vector.");
  flm_real = mxGetPr(prhs[iin]);
  flm_imag = mxGetPi(prhs[iin]);
  flm = (complex double*)malloc(flm_m * flm_n * sizeof(complex double));
  for (i=0; i<flm_m*flm_n; i++)
    flm[i] = flm_real[i];
  if( mxIsComplex(prhs[iin]) ) {
    flm_imag = mxGetPi(prhs[iin]);
    for (i=0; i<flm_m*flm_n; i++)
      flm[i] += I * flm_imag[i];
  }

  /* Parse harmonic band-limit L. */
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:bandLimit",
		      "Harmonic band-limit must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit must be positive integer.");
  if (flm_m * flm_n != L * L)
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:flmSize",
		      "Invalid number of harmonic coefficients.");

  /* Parse method. */
  iin = 2;
  if( !mxIsChar(prhs[iin]) ) {
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:methodChar",
		      "Method must be string.");
  }
  len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
  if (len >= SSHT_STRING_LEN) 
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:methodTooLong",
		      "Method exceeds string length.");
  mxGetString(prhs[iin], method, len);

  /* Parse spin. */
  iin = 3;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:spin",
		      "Spin number must be integer.");
  }
  spin = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)spin)
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:spinNonInt",
		      "Spin number must be integer.");
  if (spin >= L)
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:spinInvalid",
		      "Spin number must be strictly less than band-limit.");
  if (spin != 0 && reality)
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:spinReality",
		      "A spin function on the sphere cannot be real."); 

  /* Parse South pole flag. */
  iin = 5;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:southPoleExists",
		      "South pole flag must be logical.");
  south_pole_exists = mxIsLogicalScalarTrue(prhs[iin]);

  /* Parse North pole flag. */
  iin = 6;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:northPoleExists",
		      "North pole flag must be logical.");
  north_pole_exists = mxIsLogicalScalarTrue(prhs[iin]);

  /* Check polar interface usage valid. */
  if (south_pole_exists || north_pole_exists) {

    if (strcmp(method, SSHT_SAMPLING_MW) != 0 &&
	strcmp(method, SSHT_SAMPLING_MWSS) != 0) 
      mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:poleWithInvalidMethod",
			"Polar interfaces only supported by MW or MWSS methods.");

    if (!south_pole_exists)
      mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:noSouthPole",
			"South pole sample must be specified if North pole is.");

    if (south_pole_exists && !north_pole_exists)
      if (strcmp(method, SSHT_SAMPLING_MW) != 0)
	mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:noNorthPole",
			  "North pole must be specified for MWSS method.");

    if (south_pole_exists && north_pole_exists)
      if (strcmp(method, SSHT_SAMPLING_MWSS) != 0)
	mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:northPoleWithInvalidMethod",
			  "North pole must not be specified for MW method.");

  }

  /* Compute adjoint forward transform. */
  if (strcmp(method, SSHT_SAMPLING_MW) == 0) {

    ntheta = ssht_sampling_mw_ntheta(L);
    nphi = ssht_sampling_mw_nphi(L);

    if (south_pole_exists) {
      ntheta = ntheta - 1;
      if (reality) {
	fr = (double*)calloc(ntheta*nphi, sizeof(double));
	ssht_adjoint_mw_forward_sov_sym_real_pole(fr, &fr_sp, 
						  flm, L, 
						  dl_method, verbosity);
      }
      else {
	f = (complex double*)calloc(ntheta*nphi, sizeof(complex double));
	ssht_adjoint_mw_forward_sov_sym_pole(f, &f_sp, &phi_sp,
					     flm, L, spin, 
					     dl_method, verbosity);
      }
    }
    else {
      if (reality) {
	fr = (double*)calloc(ntheta*nphi, sizeof(double));
	ssht_adjoint_mw_forward_sov_sym_real(fr, flm, L, 
					     dl_method, verbosity);
      }
      else {
	f = (complex double*)calloc(ntheta*nphi, sizeof(complex double));
	ssht_adjoint_mw_forward_sov_sym(f, flm, L, spin, 
					dl_method, verbosity);
      }
    }

    /* mexPrintf("flm_m = %d; flm_n = %d\n", flm_m, flm_n); */
    /* mexPrintf("ntheta = %d; nphi = %d\n", ntheta, nphi); */
    /* mexPrintf("L = %d; spin = %d; reality = %d; verbosity = %d\n",  */
    /* 	      L, spin, reality, verbosity); */

  }
  else if (strcmp(method, SSHT_SAMPLING_MWSS) == 0) {
      
    ntheta = ssht_sampling_mw_ss_ntheta(L);
    nphi = ssht_sampling_mw_ss_nphi(L);

    if (south_pole_exists) {
      ntheta = ntheta - 2;
      if (reality) {
	fr = (double*)calloc(ntheta*nphi, sizeof(double));
	ssht_adjoint_mw_forward_sov_sym_ss_real_pole(fr, 
						     &fr_np, &fr_sp,
						     flm, L, 
						     dl_method, verbosity);
      }
      else {
	f = (complex double*)calloc(ntheta*nphi, sizeof(complex double));
	ssht_adjoint_mw_forward_sov_sym_ss_pole(f, &f_np, &phi_np, &f_sp, &phi_sp,
					       flm, L, spin, 
					       dl_method, verbosity);
      }
    }
    else {
      if (reality) {
	fr = (double*)calloc(ntheta*nphi, sizeof(double));
	ssht_adjoint_mw_forward_sov_sym_ss_real(fr, flm, L,
					     dl_method, verbosity);
      }
      else {
	f = (complex double*)calloc(ntheta*nphi, sizeof(complex double));
	ssht_adjoint_mw_forward_sov_sym_ss(f, flm, L, spin, 
					dl_method, verbosity);
      }
    }

  }
  else {

    mexErrMsgIdAndTxt("ssht_forward_adjoint_mex:InvalidInput:method",
		      "Method invalid.");

  } 

  /* Copy f to output. */
  if (reality) {

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(ntheta, nphi, mxREAL);
    f_real = mxGetPr(plhs[iout]);
    for(t=0; t<ntheta; t++)
      for(p=0; p<nphi; p++)
	f_real[p*ntheta + t] = fr[t*nphi + p];

    iout = 1;
    plhs[iout] =  mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[iout]) = fr_sp;

    iout = 2;
    plhs[iout] =  mxCreateDoubleMatrix(1, 1, mxREAL); 
    *mxGetPr(plhs[iout]) = phi_sp;

    iout = 3;
    plhs[iout] =  mxCreateDoubleMatrix(1, 1, mxREAL);
    *mxGetPr(plhs[iout]) = fr_np;

    iout = 4;
    plhs[iout] =  mxCreateDoubleMatrix(1, 1, mxREAL); 
    *mxGetPr(plhs[iout]) = phi_np;

  }
  else {

    iout = 0;
    plhs[iout] = mxCreateDoubleMatrix(ntheta, nphi, mxCOMPLEX);
    f_real = mxGetPr(plhs[iout]);
    f_imag = mxGetPi(plhs[iout]);    
    for(t=0; t<ntheta; t++) {
      for(p=0; p<nphi; p++) {
	f_real[p*ntheta + t] = creal(f[t*nphi + p]);
	f_imag[p*ntheta + t] = cimag(f[t*nphi + p]);
      }
    }

    iout = 1;
    plhs[iout] =  mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    *mxGetPr(plhs[iout]) = creal(f_sp);
    *mxGetPi(plhs[iout]) = cimag(f_sp);

    iout = 2;
    plhs[iout] =  mxCreateDoubleMatrix(1, 1, mxREAL); 
    *mxGetPr(plhs[iout]) = phi_sp;

    iout = 3;
    plhs[iout] =  mxCreateDoubleMatrix(1, 1, mxCOMPLEX);
    *mxGetPr(plhs[iout]) = creal(f_np);
    *mxGetPi(plhs[iout]) = cimag(f_np);

    iout = 4;
    plhs[iout] =  mxCreateDoubleMatrix(1, 1, mxREAL); 
    *mxGetPr(plhs[iout]) = phi_np;

  }

  /* Free memory. */
  free(flm);
  if (reality)
    free(fr);
  else
    free(f);

}
