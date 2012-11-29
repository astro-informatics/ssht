#ifndef SSHT_SHARP_UTILS
#define SSHT_SHARP_UTILS

#include "c_utils.h"
#include "sharp.h"
#include "sharp_almhelpers.h"
#include "sharp_geomhelpers.h"

void ssht_flm2alm_r (const complex double *flm, int L, double complex ***alm,
  sharp_alm_info **ainfo);
void ssht_flm2alm_c (const complex double *flm, int L, int spin,
  double complex ***alm, sharp_alm_info **ainfo);

void ssht_alm2flm_r (complex double *flm, int L, double complex **alm,
  sharp_alm_info *ainfo);
void ssht_alm2flm_c (complex double *flm, int L, int spin, double complex **alm,
  sharp_alm_info *ainfo);

void ssht_sharp_mw_forward_complex(complex double *flm, const complex double *f, int L, int spin);
void ssht_sharp_mw_forward_real(complex double *flm, const double *f, int L);
void ssht_sharp_mws_forward_complex(complex double *flm, const complex double *f, int L, int spin);
void ssht_sharp_mws_forward_real(complex double *flm, const double *f, int L);

#endif
