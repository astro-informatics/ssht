

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

  int i, L, spin, verbosity, f_m, f_n;
  double *flm_real, *flm_imag, *f_real, *f_imag;
  complex double *flm, *f;




  double theta, phi;
  int ntheta, nphi, n, t, p;
  double *thetas, *phis, *weights_unused;
  int len, iin = 0, iout = 0;
  char method[SSHT_STRING_LEN];

  /* Check number of arguments. */
  if(nrhs!=5) {
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:nrhs",
		      "Require five inputs.");
  }
  if(nlhs!=1) {
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidOutput:nlhs",
		      "Require one output.");
  }

  /* Parse function samples f. */
  if( !mxIsDouble(prhs[iin]) || !mxIsComplex(prhs[iin]) ) {
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:f",
		      "Function values must be complex doubles.");
  }
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);  
  f_real = mxGetPr(prhs[iin]);
  f_imag = mxGetPi(prhs[iin++]);
  f = (complex double*)malloc(f_m * f_n * sizeof(complex double));
  for(t=0; t<f_m; t++)
    for(p=0; p<f_n; p++)
      f[t*f_n + p] = f_real[p*f_m + t] + I * f_imag[p*f_m + t];

  /* Parse harmonic band-limit L. */
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:bandLimit",
		      "Harmonic band-limit must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin++]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit must be positive integer.");

  /* Parse method. */
  if( !mxIsChar(prhs[iin]) ) {
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:methodChar",
		      "Method must be string.");
  }
  len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
  if (len >= SSHT_STRING_LEN) 
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:methodTooLong",
		      "Method exceeds string length.");
  mxGetString(prhs[iin++], method, len);

  /* Parse spin. */
  spin = (int)mxGetScalar(prhs[iin++]);

  /* Parse verbosity. */
  verbosity = (int)mxGetScalar(prhs[iin++]);




  if (strcmp(method, SSHT_SAMPLING_MW) == 0) {


    ntheta = ssht_sampling_mw_ntheta(L);
    nphi = ssht_sampling_mw_nphi(L);

    if (ntheta != f_m || nphi != f_n)
      mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:inconsistentSizes",
		      "Number of function samples inconsistent with method and band-limit.");


    mexPrintf("f_m = %d; f_n = %d\n", f_m, f_n);
    mexPrintf("ntheta = %d; nphi = %d\n", ntheta, nphi);
    mexPrintf("L = %d; spin = %d; verbosity = %d\n", L, spin, verbosity);


    flm = (complex double*)calloc(L*L, sizeof(complex double));

    //    ssht_core_mw_inverse_sov_sym(f, flm, L, spin, verbosity);
    ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, verbosity);


    plhs[iout] = mxCreateDoubleMatrix(L*L, 1, mxCOMPLEX);
    flm_real = mxGetPr(plhs[iout]);
    flm_imag = mxGetPi(plhs[iout]);
    for (i=0; i<L*L; i++) {
      flm_real[i] = creal(flm[i]);
      flm_imag[i] = cimag(flm[i]);;
    }




  }
  else if (strcmp(method, SSHT_SAMPLING_MWSS) == 0) {
      
   

  }
  else if (strcmp(method, SSHT_SAMPLING_GL) == 0) {

    

  }
  else if (strcmp(method, SSHT_SAMPLING_DH) == 0) {

   
  }
  else {

    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:method",
		      "Method invalid.");

  } 

  /* Free memory. */
  /* free(flm); */
  /* free(f); */

}
