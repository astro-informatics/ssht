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

void mw2dw_real(const double *f, int L, int s, double *f_dh);
void mws2dw_real(const double *f, int L, int s, double *f_dh);

#endif
