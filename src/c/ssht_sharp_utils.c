#include "ssht_types.h"
#include "walltime_c.h"
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
  double spinsign = (spin==0) ? 1. : -1.;
  for (m=0; m<L; ++m)
    {
    double mfac=(m&1) ? -1.:1.;
    for (l=m; l<L; ++l)
      {
      complex double v1=flm[l*l+l+m], v2=mfac*flm[l*l+l-m];
      (*alm)[0][sharp_alm_index(*ainfo,l,m)]=spinsign*0.5*(v1+conj(v2));
      (*alm)[1][sharp_alm_index(*ainfo,l,m)]=-spinsign*0.5*_Complex_I*(v1-conj(v2));
      }
    }
  }

void ssht_alm2flm_c (complex double *flm, int L, int spin, double complex **alm,
  sharp_alm_info *ainfo)
  {
  int l,m;
  double spinsign = (spin==0) ? 1. : -1.;
  for (m=0; m<L; ++m)
    {
    double mfac=(m&1) ? -1:1;
    for (l=m; l<L; ++l)
      {
      complex double v1=alm[0][sharp_alm_index(ainfo,l,m)],
                    v2=alm[1][sharp_alm_index(ainfo,l,m)];
      flm[l*l+l+m]=spinsign*(v1+_Complex_I*v2);
      flm[l*l+l-m]=spinsign*(mfac*(conj(v1)+_Complex_I*conj(v2)));
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

void ssht_sharp_mw_forward_complex(complex double *flm, const complex double *f,
  int L, int spin)
  {
double t0=wallTime();
  int m,ith;
  int nphi=2*L-1;
  int nm=nphi;
  int nth_mw=L;
  int nth_mwfull=2*L-1;
  int nth_hw=2*L;
  int nth_hwfull=2*nth_hw-2;
  complex double **tmp1;
  ALLOC2D(tmp1,complex double,nth_hw,nphi);

  // FFT in phi direction
  {
  complex double *ttt=RALLOC(complex double,nphi);
  fftw_plan plan = fftw_plan_dft_1d(nphi,ttt,tmp1[0],
    FFTW_FORWARD,FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
  for (ith=0; ith<nth_mw; ++ith)
    fftw_execute_dft(plan,f+ith*nphi,tmp1[ith]);
  fftw_destroy_plan(plan);
  }

double t1=wallTime();
double tnorm=0;
  {
  double norm=1./(nth_mwfull*nphi);
  double dtheta=-SSHT_PI/nth_mwfull;
  complex double *fact=RALLOC(complex double,nth_mw);
  for (ith=0; ith<nth_mw; ++ith)
    fact[ith] = norm*cexp(_Complex_I*ith*dtheta);
  complex double *tmp=RALLOC(complex double,nth_hwfull);
  fftw_plan plan1 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_FORWARD,FFTW_ESTIMATE);
  fftw_plan plan2 = fftw_plan_dft_1d(nth_hwfull,tmp,tmp,FFTW_BACKWARD,FFTW_ESTIMATE);
  // loop over all m
  for (m=0; m<nm; ++m)
    {
    for (ith=0; ith<nth_mw; ++ith)
      tmp[ith]= tmp1[ith][m];

    // theta extension
int mreal= (m<L) ? m : m-(2*L-1);
mreal+=spin;
    int sign = (mreal&1) ? -1. : 1.;
    for (ith=nth_mw; ith<nth_mwfull; ++ith)
      tmp[ith]= sign*tmp[nth_mwfull-1-ith];

double tx=wallTime();
    fftw_execute(plan1);
    tmp[0]*=norm;
    for (ith=1; ith<nth_mw; ++ith)
      {
      tmp[ith] *= fact[ith];
      tmp[nth_mwfull-ith] *= conj(fact[ith]);
      }

    // zero padding
    for (ith=nth_mw; ith<nth_mwfull; ++ith)
      tmp[ith+nth_hwfull-nth_mwfull]=tmp[ith];
    for (ith=nth_mw; ith<nth_hwfull-nth_mw+1; ++ith)
      tmp[ith]=0.;
    fftw_execute(plan2);
tnorm+=wallTime()-tx;

    for (ith=0; ith<nth_hw; ++ith)
      tmp1[ith][m]=tmp[ith];
    }
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(tmp);
  DEALLOC(fact);
  }
double t2=wallTime();

  // FFT in phi direction
  {
  complex double *ttt=RALLOC(complex double,nphi);
  fftw_plan plan = fftw_plan_dft_1d(nphi,ttt,ttt,
    FFTW_BACKWARD,FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
  for (ith=0; ith<nth_hw; ++ith)
    fftw_execute_dft(plan,tmp1[ith],tmp1[ith]);
  fftw_destroy_plan(plan);
  }
double t3=wallTime();
printf("time overhead: %e %e %e %e\n",t1-t0,t2-t1,t3-t2,tnorm);
  sharp_geom_info *tinfo;
  sharp_make_hw_geom_info (nth_hw, nphi, 0., 2, 2*nphi, &tinfo);
  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(L-1,L-1,1,&alms);
  dcmplx **alm;
  double * fr=(double *)(tmp1[0]);
  double *frp[2];
  frp[0]=fr;
  frp[1]=fr+1;
  ALLOC2D(alm,dcmplx,2,(L*(L+1))/2);
  sharp_execute(SHARP_MAP2ALM,spin,alm,&frp[0],tinfo,alms,(spin==0)?2:1,SHARP_DP,NULL,NULL);
  ssht_alm2flm_c(flm,L,spin,alm,alms);
  DEALLOC2D(alm);
  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  DEALLOC2D(tmp1);
  }

void ssht_sharp_mw_forward_real(complex double *flm, const double *f, int L)
  {
  int m,ith;
  int nphi=2*L-1;
  int nm=L;
  int nth_mw=L;
  int nth_mwfull=2*L-1;
  int nth_hw=2*L;
  int nth_hwfull=2*nth_hw-2;
  double **tmp1;
  ALLOC2D(tmp1,double,nth_hw,2*nm);

  // FFT in phi direction
  {
  double *ttt=RALLOC(double,nphi);
  fftw_plan plan = fftw_plan_dft_r2c_1d(nphi,ttt,(complex double *)tmp1[0],
    FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
  for (ith=0; ith<nth_mw; ++ith)
    fftw_execute_dft_r2c(plan,f+ith*nphi,(complex double *)tmp1[ith]);
  fftw_destroy_plan(plan);
  }

  {
  double norm=1./(nth_mwfull*nphi);
  double dtheta=-SSHT_PI/nth_mwfull;
  complex double *fact=RALLOC(complex double,nth_mw);
  for (ith=0; ith<nth_mw; ++ith)
    fact[ith] = norm*cexp(_Complex_I*ith*dtheta);
  complex double *tmp=RALLOC(complex double,nth_hwfull);
  fftw_plan plan1 = fftw_plan_dft_1d(nth_mwfull,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan plan2 = fftw_plan_dft_1d(nth_hwfull,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
  // loop over all m
  for (m=0; m<nm; ++m)
    {
    for (ith=0; ith<nth_mw; ++ith)
      tmp[ith]= tmp1[ith][2*m]+_Complex_I*tmp1[ith][2*m+1];

    // theta extension
    int sign = (m&1) ? -1. : 1.;
    for (ith=nth_mw; ith<nth_mwfull; ++ith)
      tmp[ith]= sign*tmp[nth_mwfull-1-ith];

    fftw_execute(plan1);

    tmp[0]*=norm;
    for (ith=1; ith<nth_mw; ++ith)
      {
      tmp[ith] *= fact[ith];
      tmp[nth_mwfull-ith] *= conj(fact[ith]);
      }

    // zero padding
    for (ith=nth_mw; ith<nth_mwfull; ++ith)
      tmp[ith+nth_hwfull-nth_mwfull]=tmp[ith];
    for (ith=nth_mw; ith<nth_hwfull-nth_mw+1; ++ith)
      tmp[ith]=0.;

    fftw_execute(plan2);

    for (ith=0; ith<nth_hw; ++ith)
      {
      tmp1[ith][2*m]=creal(tmp[ith]);
      tmp1[ith][2*m+1]=cimag(tmp[ith]);
      }
    }
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(tmp);
  DEALLOC(fact);
  }

  // FFT in phi direction
  {
  complex double *ttt=RALLOC(complex double,nm);
  fftw_plan plan = fftw_plan_dft_c2r_1d(nphi,ttt,(double *)ttt,
    FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
  for (ith=0; ith<nth_hw; ++ith)
    fftw_execute_dft_c2r(plan,(complex double *)tmp1[ith],tmp1[ith]);
  fftw_destroy_plan(plan);
  }

  sharp_geom_info *tinfo;
  sharp_make_hw_geom_info (nth_hw, nphi, 0., 1, 2*nm, &tinfo);
  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(L-1,L-1,1,&alms);
  dcmplx **alm;
  ALLOC2D(alm,dcmplx,1,(L*(L+1))/2);
  sharp_execute(SHARP_MAP2ALM,0,alm,tmp1,tinfo,alms,1,SHARP_DP,NULL,NULL);
  ssht_alm2flm_r(flm,L,alm,alms);
  DEALLOC2D(alm);
  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  DEALLOC2D(tmp1);
  }

void ssht_sharp_mws_forward_complex(complex double *flm, const complex double *f, int L, int spin)
  {
  int m,ith;
  int nphi=2*L;
  int nth_hw=2*L-1;
  complex double **tmp1;
  ALLOC2D(tmp1,complex double,nth_hw,2*L);

  // FFT in phi direction
  {
  complex double *ttt=RALLOC(complex double,nphi);
  fftw_plan plan = fftw_plan_dft_1d(nphi,ttt,tmp1[0],
    FFTW_FORWARD,FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
  for (ith=0; ith<L+1; ++ith)
    fftw_execute_dft(plan,f+ith*nphi,tmp1[ith]);
  fftw_destroy_plan(plan);
  }

  {
  double norm=1./(2.*L*nphi);
  complex double *tmp=RALLOC(complex double,4*L-4);
  fftw_plan plan1 = fftw_plan_dft_1d(2*L,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan plan2 = fftw_plan_dft_1d(4*L-4,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
  // loop over all m
  for (m=0; m<nphi; ++m)
    {
    for (ith=0; ith<L+1; ++ith)
      tmp[ith]= tmp1[ith][m];

    // theta extension
int mreal= (m<=L) ? m : m-(2*L);
mreal+=spin;
    int sign = (mreal&1) ? -1. : 1.;
    for (ith=L+1; ith<2*L; ++ith)
      tmp[ith]= sign*tmp[2*L-ith];

    fftw_execute(plan1);

    for (ith=0; ith<2*L; ++ith)
      tmp[ith] *= norm;

    // zero padding
    for (ith=L+1; ith<2*L; ++ith)
      tmp[ith+2*L-4]=tmp[ith];
    for (ith=L+1; ith<3*L-3; ++ith)
      tmp[ith]=0.;

    fftw_execute(plan2);

    for (ith=0; ith<2*L-1; ++ith)
      tmp1[ith][m]=tmp[ith];
    }
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(tmp);
  }

  // FFT in phi direction
  {
  complex double *ttt=RALLOC(complex double,nphi);
  fftw_plan plan = fftw_plan_dft_1d(nphi,ttt,ttt,
    FFTW_BACKWARD,FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
  for (ith=0; ith<nth_hw; ++ith)
    fftw_execute_dft(plan,tmp1[ith],tmp1[ith]);
  fftw_destroy_plan(plan);
  }

  sharp_geom_info *tinfo;
  sharp_make_hw_geom_info (2*L-1, 2*L, 0., 2, 4*L, &tinfo);
  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(L-1,L-1,1,&alms);
  double * fr=(double *)(tmp1[0]);
  double *frp[2];
  frp[0]=fr;
  frp[1]=fr+1;
  dcmplx **alm;
  ALLOC2D(alm,dcmplx,2,(L*(L+1))/2);
  sharp_execute(SHARP_MAP2ALM,spin,alm,&frp[0],tinfo,alms,(spin==0)?2:1,SHARP_DP,NULL,NULL);
  ssht_alm2flm_c(flm,L,spin,alm,alms);
  DEALLOC2D(alm);
  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  DEALLOC2D(tmp1);
  }

void ssht_sharp_mws_forward_real(complex double *flm, const double *f, int L)
  {
  int m,ith;
  int nphi=2*L;
  int nth_hw=2*L-1;
  double **tmp1;
  ALLOC2D(tmp1,double,nth_hw,2*L+2);

  // FFT in phi direction
  {
  double *ttt=RALLOC(double,nphi);
  fftw_plan plan = fftw_plan_dft_r2c_1d(nphi,ttt,(complex double *)tmp1[0],
    FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
  for (ith=0; ith<L+1; ++ith)
    fftw_execute_dft_r2c(plan,f+ith*nphi,(complex double *)tmp1[ith]);
  fftw_destroy_plan(plan);
  }

  {
  double norm=1./(2.*L*nphi);
  complex double *tmp=RALLOC(complex double,4*L-4);
  fftw_plan plan1 = fftw_plan_dft_1d(2*L,tmp,tmp,FFTW_FORWARD,FFTW_MEASURE);
  fftw_plan plan2 = fftw_plan_dft_1d(4*L-4,tmp,tmp,FFTW_BACKWARD,FFTW_MEASURE);
  // loop over all m
  for (m=0; m<L+1; ++m)
    {
    for (ith=0; ith<L+1; ++ith)
      tmp[ith]= tmp1[ith][2*m]+_Complex_I*tmp1[ith][2*m+1];

    // theta extension
    int sign = (m&1) ? -1. : 1.;
    for (ith=L+1; ith<2*L; ++ith)
      tmp[ith]= sign*tmp[2*L-ith];

    fftw_execute(plan1);

    for (ith=0; ith<2*L; ++ith)
      tmp[ith] *= norm;

    // zero padding
    for (ith=L+1; ith<2*L; ++ith)
      tmp[ith+2*L-4]=tmp[ith];
    for (ith=L+1; ith<3*L-3; ++ith)
      tmp[ith]=0.;

    fftw_execute(plan2);

    for (ith=0; ith<2*L-1; ++ith)
      {
      tmp1[ith][2*m]=creal(tmp[ith]);
      tmp1[ith][2*m+1]=cimag(tmp[ith]);
      }
    }
  fftw_destroy_plan(plan1);
  fftw_destroy_plan(plan2);
  DEALLOC(tmp);
  }

  // FFT in phi direction
  {
  complex double *ttt=RALLOC(complex double,L+1);
  fftw_plan plan = fftw_plan_dft_c2r_1d(nphi,ttt,(double *)ttt,
    FFTW_MEASURE|FFTW_UNALIGNED);
  DEALLOC(ttt);
  for (ith=0; ith<nth_hw; ++ith)
    fftw_execute_dft_c2r(plan,(complex double *)tmp1[ith],tmp1[ith]);
  fftw_destroy_plan(plan);
  }

  sharp_geom_info *tinfo;
  sharp_make_hw_geom_info (2*L-1, 2*L, 0., 1, 2*L+2, &tinfo);
  sharp_alm_info *alms;
  sharp_make_triangular_alm_info(L-1,L-1,1,&alms);
  dcmplx **alm;
  ALLOC2D(alm,dcmplx,1,(L*(L+1))/2);
  sharp_execute(SHARP_MAP2ALM,0,alm,tmp1,tinfo,alms,1,SHARP_DP,NULL,NULL);
  ssht_alm2flm_r(flm,L,alm,alms);
  DEALLOC2D(alm);
  sharp_destroy_alm_info(alms);
  sharp_destroy_geom_info(tinfo);

  DEALLOC2D(tmp1);
  }
