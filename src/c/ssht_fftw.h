// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SSHT_FFTW
#define SSHT_FFTW

#ifdef __cplusplus
#include "ssht_types.h"
#include <fftw3.h>
fftw_plan fftw_plan_dft_2d(int n0, int n1, std::complex<double> * input, std::complex<double> * output,
			int direction, unsigned int method);
fftw_plan fftw_plan_dft_c2r_2d(int n0, int n1, std::complex<double> * input, double * output,
			 unsigned int method);
fftw_plan fftw_plan_dft_c2r_1d(int n0, std::complex<double> * input, double * output,
			 unsigned int method);
fftw_plan fftw_plan_dft_1d(int n0, std::complex<double> * input, std::complex<double> * output, int direction, unsigned int method);
fftw_plan fftw_plan_dft_r2c_2d(int n0, int n1, double * input,std::complex<double> * output,
			unsigned int method);
fftw_plan fftw_plan_dft_r2c_1d(int n0, double * input, std::complex<double> * output,
			unsigned int method);
//! exectute functions for FFTW
void fftw_execute_dft(fftw_plan plan, std::complex<double> * in, std::complex<double> * out);
void fftw_execute_dft_c2r(fftw_plan plan, std::complex<double> * in, double * out);
void fftw_execute_dft_r2c(fftw_plan plan, double * in, std::complex<double> * out);
#endif
#endif
