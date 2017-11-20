// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details



#include "ssht_fftw.h"
#ifdef __cplusplus
fftw_plan fftw_plan_dft_2d(int n0, int n1, std::complex<double> * input, std::complex<double> * output,
			int direction, unsigned int method){
    return fftw_plan_dft_2d(n0, n1, reinterpret_cast<fftw_complex *>(input), reinterpret_cast<fftw_complex *>(output),
        direction, method);
}
fftw_plan fftw_plan_dft_c2r_2d(int n0, int n1, std::complex<double> * input, double * output,
			 unsigned int method){
return fftw_plan_dft_c2r_2d(n0, n1, reinterpret_cast<fftw_complex *>(input), output,
			      method);
}
fftw_plan fftw_plan_dft_c2r_1d(int n0, std::complex<double> * input, double * output,
			 unsigned int method){
return fftw_plan_dft_c2r_1d(n0, reinterpret_cast<fftw_complex *>(input), output,
			      method);
}
fftw_plan fftw_plan_dft_1d(int n0, std::complex<double> * input, std::complex<double> * output, int direction, unsigned int method)
{
  return fftw_plan_dft_1d(n0, reinterpret_cast<fftw_complex *>(input), reinterpret_cast<fftw_complex *>(output), direction, method);

}
fftw_plan fftw_plan_dft_r2c_2d(int n0, int n1, double * input,std::complex<double> * output,
			unsigned int method){
return fftw_plan_dft_r2c_2d(n0, n1,input, reinterpret_cast<fftw_complex *>(output),
			      method);
}
fftw_plan fftw_plan_dft_r2c_1d(int n0, double * input, std::complex<double> * output,
			unsigned int method){
return fftw_plan_dft_r2c_1d(n0,input, reinterpret_cast<fftw_complex *>(output),
			      method);
}
//! exectute functions for FFTW
void fftw_execute_dft(fftw_plan plan, std::complex<double> * in, std::complex<double> * out){
fftw_execute_dft(plan,reinterpret_cast<std::complex<double> *>(in), reinterpret_cast<std::complex<double> *>(out));
}
void fftw_execute_dft_c2r(fftw_plan plan, std::complex<double> * in, double * out){
fftw_execute_dft_c2r(plan, reinterpret_cast<std::complex<double> *>(in), out);
}
void fftw_execute_dft_r2c(fftw_plan plan, double * in, std::complex<double> * out){
fftw_execute_dft_r2c(plan, in, reinterpret_cast<std::complex<double> *>(out));
}
#endif
