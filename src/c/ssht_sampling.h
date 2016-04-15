// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


#ifndef SSHT_SAMPLING
#define SSHT_SAMPLING

#include <complex.h>

complex double ssht_sampling_weight_mw(int p);
double ssht_sampling_weight_dh(double theta_t, int L);
void ssht_sampling_gl_thetas_weights(double *thetas, double *weights, int L);

double ssht_sampling_mw_t2theta(int t, int L);
double ssht_sampling_mw_p2phi(int p, int L);
int ssht_sampling_mw_n(int L);
int ssht_sampling_mw_ntheta(int L);
int ssht_sampling_mw_nphi(int L);

double ssht_sampling_mw_ss_t2theta(int t, int L);
double ssht_sampling_mw_ss_p2phi(int p, int L);
int ssht_sampling_mw_ss_n(int L);
int ssht_sampling_mw_ss_ntheta(int L);
int ssht_sampling_mw_ss_nphi(int L);

double ssht_sampling_dh_t2theta(int t, int L);
double ssht_sampling_dh_p2phi(int p, int L);
int ssht_sampling_dh_n(int L);
int ssht_sampling_dh_ntheta(int L);
int ssht_sampling_dh_nphi(int L);

double ssht_sampling_gl_p2phi(int p, int L);
int ssht_sampling_gl_n(int L);
int ssht_sampling_gl_ntheta(int L);
int ssht_sampling_gl_nphi(int L);

void ssht_sampling_elm2ind(int *ind, int el, int m);
void ssht_sampling_ind2elm(int *el, int *m, int ind);

#endif
