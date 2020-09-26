#include <complex.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>

#include "ssht.h"
#include "utilities.h"

#include <cmocka.h>

struct Parameters {
  int L, L0, seed, verbosity, dl_method;
};
typedef void (*real_inverse_transform)(
    double *f, const ssht_complex_double *flm, const struct Parameters *const params);
typedef void (*real_forward_transform)(
    ssht_complex_double *flm, const double *f, const struct Parameters *const params);

static int inverse_coeff_size(int L) {
  int value = L * L;
  if (L * (2 * L - 1) > value)
    value = L * (2 * L - 1);
  return value;
}

static void back_and_forth(
    real_inverse_transform inverse,
    real_forward_transform forward,
    const struct Parameters *const params) {
  complex double *complex_original =
      (complex double *)calloc(inverse_coeff_size(params->L), sizeof(complex double));
  double *real =
      (double *)calloc((2 * params->L) * (2 * params->L - 1), sizeof(double));
  complex double *complex_final =
      (complex double *)calloc(params->L * params->L, sizeof(complex double));

  (*inverse)(real, complex_original, params);
  (*forward)(complex_final, real, params);

  for (int i = 0; i < params->L * params->L; i++) {
    assert_float_equal((float)complex_final[i], (float)complex_original[i], 1e-12);
  }

  free(complex_original);
  free(real);
  free(complex_final);
}

static void dh_inverse_sov_real(
    double *f, const ssht_complex_double *flm, const struct Parameters *const params) {
  ssht_core_dh_inverse_sov_real(f, flm, params->L, params->verbosity);
}
static void dh_forward_sov_real(
    ssht_complex_double *flm, const double *f, const struct Parameters *const params) {
  ssht_core_dh_forward_sov_real(flm, f, params->L, params->verbosity);
}

static void back_and_forth_dh_real(void **state) {
  back_and_forth(
      dh_inverse_sov_real,
      dh_forward_sov_real,
      (const struct Parameters *const) * state);
}

static void gl_inverse_sov_real(
    double *f, const ssht_complex_double *flm, const struct Parameters *const params) {
  ssht_core_gl_inverse_sov_real(f, flm, params->L, params->verbosity);
}
static void gl_forward_sov_real(
    ssht_complex_double *flm, const double *f, const struct Parameters *const params) {
  ssht_core_gl_forward_sov_real(flm, f, params->L, params->verbosity);
}
static void back_and_forth_gl_real(void **state) {
  back_and_forth(
      gl_inverse_sov_real,
      gl_forward_sov_real,
      (const struct Parameters *const) * state);
}

static void sym_inverse_sov_real(
    double *f, const ssht_complex_double *flm, const struct Parameters *const params) {
  ssht_core_mw_inverse_sov_sym_real(
      f, flm, params->L, params->dl_method, params->verbosity);
}
static void sym_forward_sov_conv_real(
    ssht_complex_double *flm, const double *f, const struct Parameters *const params) {
  ssht_core_mw_forward_sov_conv_sym_real(
      flm, f, params->L, params->dl_method, params->verbosity);
}

static void back_and_forth_sym_real(void **state) {
  back_and_forth(
      sym_inverse_sov_real,
      sym_forward_sov_conv_real,
      (const struct Parameters *const) * state);
}

static void mw_lb_inverse_sov_sym_real(
    double *f, const ssht_complex_double *flm, const struct Parameters *const params) {
  ssht_core_mw_lb_inverse_sov_sym_real(
      f, flm, params->L0, params->L, params->dl_method, params->verbosity);
}
static void mw_lb_forward_sov_conv_sym_real(
    ssht_complex_double *flm, const double *f, const struct Parameters *const params) {
  ssht_core_mw_lb_forward_sov_conv_sym_real(
      flm, f, params->L0, params->L, params->dl_method, params->verbosity);
}

static void back_and_forth_lb_real(void **state) {
  back_and_forth(
      mw_lb_inverse_sov_sym_real,
      mw_lb_forward_sov_conv_sym_real,
      (const struct Parameters *const) * state);
}

int main(void) {
  const struct Parameters params[] = {
      {.L = 64, .L0 = 32, .seed = 1, .verbosity = 0, .dl_method = SSHT_DL_RISBO},
      {.L = 64, .L0 = 32, .seed = 1, .verbosity = 0, .dl_method = SSHT_DL_TRAPANI}};
  const struct CMUnitTest tests[] = {
      cmocka_unit_test_prestate(back_and_forth_dh_real, (void *)params),
      cmocka_unit_test_prestate(back_and_forth_gl_real, (void *)params),
      cmocka_unit_test_prestate(back_and_forth_sym_real, (void *)params),
      cmocka_unit_test_prestate(back_and_forth_sym_real, (void *)(params + 1)),
      cmocka_unit_test_prestate(back_and_forth_lb_real, (void *)params),
      cmocka_unit_test_prestate(back_and_forth_lb_real, (void *)(params + 1)),
  };

  return cmocka_run_group_tests(tests, NULL, NULL);
}
