

#include <ssht.h>
#include "ssht_mex.h"
#include <string.h>
#include "mex.h"


/**
 * Compute inverse transform.
 *
 * Usage: [n, ntheta, nphi, thetas, phis] = ssht_sampling_mex(L, method)




 *
 * \author Jason McEwen
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  int i, L, spin, reality, verbosity=0, flm_m, flm_n;
  double *flm_real, *flm_imag, *f_real, *f_imag;
  complex double *flm, *f;
  int ntheta, nphi, t, p;
  int len, iin = 0, iout = 0;
  char method[SSHT_STRING_LEN];

  /* Check number of arguments. */
  if(nrhs!=5) {
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:nrhs",
		      "Require five inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidOutput:nlhs",
		      "Require one output.");
  }

  /* Parse harmonic coefficients flm. */
  if( !mxIsDouble(prhs[iin]) || !mxIsComplex(prhs[iin]) ) {
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:flm",
		      "Harmonic coefficients must be complex doubles.");
  }
  flm_m = mxGetM(prhs[iin]);
  flm_n = mxGetN(prhs[iin]);
  if (flm_m != 1 && flm_n != 1)
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:flmVector",
		      "Harmonic coefficients must be contained in vector.");
  flm_real = mxGetPr(prhs[iin]);
  flm_imag = mxGetPi(prhs[iin++]);
  flm = (complex double*)malloc(flm_m * flm_n * sizeof(complex double));
  for (i=0; i<flm_m*flm_n; i++)
    flm[i] = flm_real[i] + I * flm_imag[i];

  /* Parse harmonic band-limit L. */
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:bandLimit",
		      "Harmonic band-limit must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin++]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit must be positive integer.");
  if (flm_m * flm_n != L * L)
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:flmSize",
		      "Invalid number of harmonic coefficients.");

  /* Parse method. */
  if( !mxIsChar(prhs[iin]) ) {
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:methodChar",
		      "Method must be string.");
  }
  len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
  if (len >= SSHT_STRING_LEN) 
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:methodTooLong",
		      "Method exceeds string length.");
  mxGetString(prhs[iin++], method, len);

  /* Parse spin. */
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:spin",
		      "Spin number must be integer.");
  }
  spin = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin++]) > (double)spin || spin < 0)
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:spinNonInt",
		      "Spin number must be positive integer.");
  if (spin >= L)
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:spinInvalid",
		      "Spin number must be strictly less than band-limit.");

  /* Parse reality. */
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:reality",
		      "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin++]);

  /* Compute inverse transform. */
  if (strcmp(method, SSHT_SAMPLING_MW) == 0) {

    ntheta = ssht_sampling_mw_ntheta(L);
    nphi = ssht_sampling_mw_nphi(L);
    f = (complex double*)calloc(ntheta*nphi, sizeof(complex double));
    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, verbosity);

    /* mexPrintf("flm_m = %d; flm_n = %d\n", flm_m, flm_n); */
    /* mexPrintf("ntheta = %d; nphi = %d\n", ntheta, nphi); */
    /* mexPrintf("L = %d; spin = %d; reality = %d; verbosity = %d\n",  */
    /* 	      L, spin, reality, verbosity); */

  }
  else if (strcmp(method, SSHT_SAMPLING_MWSS) == 0) {
      
    ntheta = ssht_sampling_mw_ss_ntheta(L);
    nphi = ssht_sampling_mw_ss_nphi(L);
    f = (complex double*)calloc(ntheta*nphi, sizeof(complex double));
    ssht_core_mw_inverse_sov_sym_ss(f, flm, L, spin, verbosity);

  }
  else if (strcmp(method, SSHT_SAMPLING_GL) == 0) {

    ntheta = ssht_sampling_gl_ntheta(L);
    nphi = ssht_sampling_gl_nphi(L);
    f = (complex double*)calloc(ntheta*nphi, sizeof(complex double));
    ssht_core_gl_inverse_sov(f, flm, L, spin, verbosity);

  }
  else if (strcmp(method, SSHT_SAMPLING_DH) == 0) {

    ntheta = ssht_sampling_dh_ntheta(L);
    nphi = ssht_sampling_dh_nphi(L);
    f = (complex double*)calloc(ntheta*nphi, sizeof(complex double));
    ssht_core_dh_inverse_sov(f, flm, L, spin, verbosity);
   
  }
  else {

    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:method",
		      "Method invalid.");

  } 

  /* Copy f to output. */
  plhs[iout] = mxCreateDoubleMatrix(ntheta, nphi, mxCOMPLEX);
  f_real = mxGetPr(plhs[iout]);
  f_imag = mxGetPi(plhs[iout]);    
  for(t=0; t<ntheta; t++) {
    for(p=0; p<nphi; p++) {
      f_real[p*ntheta + t] = creal(f[t*nphi + p]);
      f_imag[p*ntheta + t] = cimag(f[t*nphi + p]);
    }
  }

  /* Free memory. */
  free(flm);
  free(f);

}
