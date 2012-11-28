#include "ssht_types.h"
#include "ssht_sharp_utils.h"
#include <fftw3.h>

typedef double complex dcmplx;

void ssht_flm2alm_r (const complex double *flm, int L, double complex ***alm,
  sharp_alm_info **ainfo)
  {
  sharp_make_triangular_alm_info(L-1,L-1,1,ainfo);
  ALLOC2D((*alm),dcmplx,1,(L*(L+1))/2);
  // rearrange a_lm into required format
  int l,m;
  for (m=0; m<L; ++m)
    for (l=m; l<L; ++l)
      (*alm)[0][sharp_alm_index(*ainfo,l,m)]=flm[l*l + l + m];
  }

void ssht_flm2alm_c (const complex double *flm, int L, int spin,
  double complex ***alm, sharp_alm_info **ainfo)
  {
  sharp_make_triangular_alm_info(L-1,L-1,1,ainfo);

  ALLOC2D((*alm),dcmplx,2,(L*(L+1))/2);
  // rearrange a_lm into required format
  int l,m;
  for (m=0; m<L; ++m)
    {
    double mfac=(m&1) ? -1.:1.;
    for (l=m; l<L; ++l)
      {
      complex double v1=flm[l*l+l+m], v2=mfac*flm[l*l+l-m];
      (*alm)[0][sharp_alm_index(*ainfo,l,m)]=0.5*(v1+conj(v2));
      (*alm)[1][sharp_alm_index(*ainfo,l,m)]=-0.5*_Complex_I*(v1-conj(v2));
      }
    }
  }

void ssht_alm2flm_c (complex double *flm, int L, int spin, double complex **alm,
  sharp_alm_info *ainfo)
  {
  int l,m;
  for (m=0; m<L; ++m)
    {
    double mfac=(m&1) ? -1:1;
    for (l=m; l<L; ++l)
      {
      complex double v1=alm[0][sharp_alm_index(ainfo,l,m)],
                     v2=alm[1][sharp_alm_index(ainfo,l,m)];
      flm[l*l+l+m]=v1+_Complex_I*v2;
      flm[l*l+l-m]=mfac*(conj(v1)+_Complex_I*conj(v2));
      }
    }
  }

void ssht_alm2flm_r (complex double *flm, int L, double complex **alm,
  sharp_alm_info *ainfo)
  {
  int l,m;
  for (m=0; m<L; ++m)
    for (l=m; l<L; ++l)
      {
      flm[l*l + l + m]=alm[0][sharp_alm_index(ainfo,l,m)];
      flm[l*l + l - m]=conj(flm[l*l + l + m])*((m&1)? -1:1);
      }
  }

void mw2dw_real(const double *f, int L, int s, double *f_dh)
  {
  int m,ith;
  int nphi=2*L-1;
  int nth_dh=2*L;
  complex double **tmp1;
  ALLOC2D(tmp1,complex double,nth_dh,L);

  // FFT in phi direction
  {
  double *ttt=RALLOC(double,nphi);
  fftw_plan plan = fftw_plan_dft_r2c_1d(nphi,ttt,tmp1[0],
    FFTW_MEASURE|FFTW_UNALIGNED);
  for (ith=0; ith<L; ++ith)
    fftw_execute_dft_r2c(plan,f+ith*nphi,tmp1[ith]);
  fftw_destroy_plan(plan);
  DEALLOC(ttt);
  }

  {
double norm=1./((2.*L-1.)*nphi);
double dtheta=SSHT_PI*(-1./(2.0*L-1.0)+1./(4.*L));
  complex double *tmp=RALLOC(complex double,4*L);
  fftw_plan plan1 = fftw_plan_dft_1d(2*L-1,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan plan2 = fftw_plan_dft_1d(4*L,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
  // loop over all m
  for (m=0; m<L; ++m)
    {
    for (ith=0; ith<L; ++ith)
      tmp[ith]= tmp1[ith][m];

    // theta extension
    int sign = ((m+s)&1) ? -1. : 1.;
    for (ith=L; ith<2*L-1; ++ith)
      tmp[ith]= sign*tmp[2*L-2-ith];

    fftw_execute(plan1);

    tmp[0]*=norm;
    for (ith=1; ith<L; ++ith)
      {
      tmp[ith] *= norm*cexp(_Complex_I*ith*dtheta);
      tmp[2*L-1-ith] *= norm*cexp(-_Complex_I*ith*dtheta);
      }

    // zero padding
    for (ith=L; ith<2*L-1; ++ith)
      tmp[ith+2*L+1]=tmp[ith];
    for (ith=L; ith<3*L+1; ++ith)
      tmp[ith]=0.;

    fftw_execute(plan2);

    for (ith=0; ith<2*L; ++ith)
      tmp1[ith][m]=tmp[ith];
    }
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(tmp);
  }

  // FFT in phi direction
  {
  complex double *ttt=RALLOC(complex double,L);
  fftw_plan plan = fftw_plan_dft_c2r_1d(nphi,ttt,f_dh,
    FFTW_MEASURE|FFTW_UNALIGNED);
  for (ith=0; ith<nth_dh; ++ith)
    fftw_execute_dft_c2r(plan,tmp1[ith],f_dh+ith*nphi);
  fftw_destroy_plan(plan);
  DEALLOC(ttt);
  }

  DEALLOC2D(tmp1);
  }

void mws2dw_real(const double *f, int L, int s, double *f_dh)
  {
  int m,ith;
  int nphi=2*L;
  int nth_dh=2*L;
  int nth_mws=L+1;
  complex double **tmp1;
  ALLOC2D(tmp1,complex double,nth_dh,L+1);

  // FFT in phi direction
  {
  double *ttt=RALLOC(double,nphi);
  fftw_plan plan = fftw_plan_dft_r2c_1d(nphi,ttt,tmp1[0],
    FFTW_MEASURE|FFTW_UNALIGNED);
  for (ith=0; ith<L+1; ++ith)
    fftw_execute_dft_r2c(plan,f+ith*nphi,tmp1[ith]);
  fftw_destroy_plan(plan);
  DEALLOC(ttt);
  }

  {
double norm=1./((2.*L)*nphi);

  complex double *tmp=RALLOC(complex double,4*L);
  fftw_plan plan1 = fftw_plan_dft_1d(2*L,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan plan2 = fftw_plan_dft_1d(4*L,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
  // loop over all m
  for (m=0; m<L; ++m)
    {
    for (ith=0; ith<nth_mws; ++ith)
      tmp[ith]= tmp1[ith][m];
    // theta extension
    int sign = ((m+s)&1) ? -1. : 1.;
    for (ith=L+1; ith<2*L; ++ith)
      tmp[ith]= sign*tmp[2*L-ith];

    fftw_execute(plan1);

    for (ith=0; ith<2*L; ++ith)
      tmp[ith] *= norm;

    // zero padding
    for (ith=L+1; ith<2*L; ++ith)
      tmp[ith+2*L]=tmp[ith];
    for (ith=L+1; ith<3*L; ++ith)
      tmp[ith]=0.;

    fftw_execute(plan2);

    for (ith=0; ith<2*L; ++ith)
      tmp1[ith][m]=tmp[ith];
    }
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(tmp);
  }

  // FFT in phi direction
  {
  complex double *ttt=RALLOC(complex double,L+1);
  fftw_plan plan = fftw_plan_dft_c2r_1d(nphi,ttt,f_dh,
    FFTW_MEASURE|FFTW_UNALIGNED);
  for (ith=0; ith<nth_dh; ++ith)
    fftw_execute_dft_c2r(plan,tmp1[ith],f_dh+ith*nphi);
  fftw_destroy_plan(plan);
  DEALLOC(ttt);
  }

  DEALLOC2D(tmp1);
  }
