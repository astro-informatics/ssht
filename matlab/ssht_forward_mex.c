

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

  int i, L, spin, reality, verbosity=0, f_m, f_n;
  int f_is_complex;
  double *flm_real, *flm_imag, *f_real, *f_imag;
  double *fr;
  complex double *flm, *f;
  int ntheta, nphi, t, p;
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

  /* Parse reality. */
  iin = 4;
  if( !mxIsLogicalScalar(prhs[iin]) )
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:reality",
		      "Reality flag must be logical.");
  reality = mxIsLogicalScalarTrue(prhs[iin]);
  if (!f_is_complex && !reality)
    mexWarnMsgTxt("Running complex transform on real signal (set reality flag to improve performance).");

  /* Parse function samples f. */
  iin = 0;
  if( !mxIsDouble(prhs[iin]) ) {
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:f",
		      "Function values must be doubles.");
  }
  f_m = mxGetM(prhs[iin]);
  f_n = mxGetN(prhs[iin]);  
  f_real = mxGetPr(prhs[iin]);  
  f_is_complex = mxIsComplex(prhs[iin]);
  f_imag = f_is_complex ? mxGetPi(prhs[iin]) : NULL;
  if (reality) {
    fr = (double*)malloc(f_m * f_n * sizeof(double));
    for(t=0; t<f_m; t++)
      for(p=0; p<f_n; p++)
	fr[t*f_n + p] = f_real[p*f_m + t];
  }
  else {
    f = (complex double*)malloc(f_m * f_n * sizeof(complex double));
    for(t=0; t<f_m; t++)
      for(p=0; p<f_n; p++)
	f[t*f_n + p] = f_real[p*f_m + t] 
	  + I * (f_is_complex ? f_imag[p*f_m + t] : 0.0);
  }
  if (f_is_complex && reality)
    mexWarnMsgTxt("Running real transform but input appears to be complex (ignoring imaginary component).");

  /* Parse harmonic band-limit L. */
  iin = 1;
  if( !mxIsDouble(prhs[iin]) || 
      mxIsComplex(prhs[iin]) || 
      mxGetNumberOfElements(prhs[iin])!=1 ) {
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:bandLimit",
		      "Harmonic band-limit must be integer.");
  }
  L = (int)mxGetScalar(prhs[iin]);
  if (mxGetScalar(prhs[iin]) > (double)L || L <= 0)
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:bandLimitNonInt",
		      "Harmonic band-limit must be positive integer.");

  /* Parse method. */
  iin = 2;
  if( !mxIsChar(prhs[iin]) ) {
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:methodChar",
		      "Method must be string.");
  }
  len = (mxGetM(prhs[iin]) * mxGetN(prhs[iin])) + 1;
  if (len >= SSHT_STRING_LEN) 
    mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:methodTooLong",
		      "Method exceeds string length.");
  mxGetString(prhs[iin], method, len);

  /* Parse spin. */
  iin = 3;
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
  if (spin != 0 && reality)
    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:spinReality",
		      "A spin function on the sphere cannot be real."); 

  /* Compute forward transform. */
  if (strcmp(method, SSHT_SAMPLING_MW) == 0) {

    ntheta = ssht_sampling_mw_ntheta(L);
    nphi = ssht_sampling_mw_nphi(L);
    if (ntheta != f_m || nphi != f_n)
      mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:inconsistentSizesMW",
        "Number of function samples inconsistent with method and band-limit.");
    flm = (complex double*)calloc(L*L, sizeof(complex double));
    if (reality)
      ssht_core_mw_forward_sov_conv_sym_real(flm, fr, L, verbosity);
    else
      ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, verbosity);

    /* mexPrintf("flm_m = %d; flm_n = %d\n", flm_m, flm_n); */
    /* mexPrintf("ntheta = %d; nphi = %d\n", ntheta, nphi); */
    /* mexPrintf("L = %d; spin = %d; reality = %d; verbosity = %d\n",  */
    /* 	      L, spin, reality, verbosity); */

  }
  else if (strcmp(method, SSHT_SAMPLING_MWSS) == 0) {
      
    ntheta = ssht_sampling_mw_ss_ntheta(L);
    nphi = ssht_sampling_mw_ss_nphi(L);
    if (ntheta != f_m || nphi != f_n)
      mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:inconsistentSizesMWSS",
        "Number of function samples inconsistent with method and band-limit.");
    flm = (complex double*)calloc(L*L, sizeof(complex double));
    if (reality)
      ssht_core_mw_forward_sov_conv_sym_ss_real(flm, fr, L, verbosity);   
    else
      ssht_core_mw_forward_sov_conv_sym_ss(flm, f, L, spin, verbosity);   

  }
  else if (strcmp(method, SSHT_SAMPLING_GL) == 0) {

    ntheta = ssht_sampling_gl_ntheta(L);
    nphi = ssht_sampling_gl_nphi(L);
    if (ntheta != f_m || nphi != f_n)
      mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:inconsistentSizesGL",
        "Number of function samples inconsistent with method and band-limit.");
    flm = (complex double*)calloc(L*L, sizeof(complex double));
    if (reality)
      ssht_core_gl_forward_sov_real(flm, fr, L, verbosity);    
    else
      ssht_core_gl_forward_sov(flm, f, L, spin, verbosity);    

  }
  else if (strcmp(method, SSHT_SAMPLING_DH) == 0) {

    ntheta = ssht_sampling_dh_ntheta(L);
    nphi = ssht_sampling_dh_nphi(L);
    if (ntheta != f_m || nphi != f_n)
      mexErrMsgIdAndTxt("ssht_forward_mex:InvalidInput:inconsistentSizesDH",
        "Number of function samples inconsistent with method and band-limit.");
    flm = (complex double*)calloc(L*L, sizeof(complex double));
    if (reality)
      ssht_core_dh_forward_sov_real(flm, fr, L, verbosity);    
    else
      ssht_core_dh_forward_sov(flm, f, L, spin, verbosity);    
   
  }
  else {

    mexErrMsgIdAndTxt("ssht_inverse_mex:InvalidInput:method",
		      "Method invalid.");

  } 

  /* Copy flm to output. */
  plhs[iout] = mxCreateDoubleMatrix(L*L, 1, mxCOMPLEX);
  flm_real = mxGetPr(plhs[iout]);
  flm_imag = mxGetPi(plhs[iout]);
  for (i=0; i<L*L; i++) {
    flm_real[i] = creal(flm[i]);
    flm_imag[i] = cimag(flm[i]);;
  }

  /* Free memory. */
  free(flm);
  if (reality)
    free(fr);
  else
    free(f);
}
