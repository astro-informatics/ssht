#include <complex.h>
#include <setjmp.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "interface.h"
#include "ssht.h"
#include "ssht_error.h"
#include "utilities.h"

#include <cmocka.h>

struct State {
  struct ssht_InterfaceParameters params;
  int seed;
};

static void real_back_and_forth(void **void_state) {
  const struct State *const state = (const struct State *const) * void_state;
  const struct ssht_InterfaceParameters *const params = &state->params;
  complex double *complex_original =
      (complex double *)calloc(harmonic_space_size(*params), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(complex_original);
  gen_flm_real_interface(complex_original, state->seed, *params);
  double *real = (double *)calloc(real_image_space_size(*params), sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(real);
  complex double *complex_final =
      (complex double *)calloc(harmonic_space_size(*params), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(complex_final);

  ssht_real_inverse(real, complex_original, *params);
  ssht_real_forward(complex_final, real, *params);

  for (int i = 0; i < params->L * params->L; i++) {
    assert_float_equal(
        (float)creal(complex_final[i]), (float)creal(complex_original[i]), 1e-12);
    assert_float_equal(
        (float)cimag(complex_final[i]), (float)cimag(complex_original[i]), 1e-12);
  }

  free(complex_original);
  free(real);
  free(complex_final);
}

static void spin_back_and_forth(void **void_state) {
  const struct State *const state = (const struct State *const) * void_state;
  const struct ssht_InterfaceParameters *const params = &state->params;
  complex double *complex_original =
      (complex double *)calloc(harmonic_space_size(*params), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(complex_original);
  gen_flm_complex_interface(complex_original, state->seed, *params);
  complex double *image =
      (complex double *)calloc(spin_image_space_size(*params), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(image);
  complex double *complex_final =
      (complex double *)calloc(harmonic_space_size(*params), sizeof(complex double));
  SSHT_ERROR_MEM_ALLOC_CHECK(complex_final);

  ssht_spin_inverse(image, complex_original, *params);
  ssht_spin_forward(complex_final, image, *params);

  for (int i = 0; i < params->L * params->L; i++) {
    assert_float_equal(
        (float)creal(complex_final[i]), (float)creal(complex_original[i]), 1e-12);
    assert_float_equal(
        (float)cimag(complex_final[i]), (float)cimag(complex_original[i]), 1e-12);
  }

  free(complex_original);
  free(image);
  free(complex_final);
}

struct CMUnitTest
parametrization(int spin, ssht_transforms method, ssht_dl_method_t dl) {
  struct State *state = (struct State *)calloc(1, sizeof(struct State));
  SSHT_ERROR_MEM_ALLOC_CHECK(state);
  state->params.L = 64;
  state->params.L0 = 32;
  state->params.spin = spin;
  state->params.verbosity = 0;
  state->params.dl_method = dl;
  state->params.method = method;
  state->seed = 1;
  char mname[256];
  method_name(mname, state->params);
  char *name = (char *)malloc(strlen(mname) + 50);
  sprintf(name, "spin % i back and forth: %s", spin, mname);
  struct CMUnitTest result = {
      .name = name,
      .test_func = spin == 0 ? real_back_and_forth : spin_back_and_forth,
      .setup_func = NULL,
      .teardown_func = NULL,
      .initial_state = state,
  };
  return result;
}

void delete_parametrization(struct CMUnitTest const *parametrization) {
  if (parametrization == NULL)
    return;
  if (parametrization->initial_state == NULL)
    free(parametrization->initial_state);
  if (parametrization->name == NULL)
    free((void *)parametrization->name);
}

int main(void) {
  const struct CMUnitTest tests[] = {
      parametrization(0, DH_SOV, SSHT_DL_RISBO),
      parametrization(0, GL_SOV, SSHT_DL_RISBO),
      parametrization(0, MW_SOV_SYM, SSHT_DL_RISBO),
      parametrization(0, MW_SOV_SYM, SSHT_DL_TRAPANI),
      parametrization(0, MW_SOV_SYM_LB, SSHT_DL_RISBO),
      parametrization(0, MW_SOV_SYM_LB, SSHT_DL_TRAPANI),
      parametrization(0, MW_SOV_SYM_POLE, SSHT_DL_RISBO),
      parametrization(0, MW_SOV_SYM_POLE, SSHT_DL_TRAPANI),
      parametrization(0, MW_SOV_SYM_LB_SS, SSHT_DL_RISBO),
      parametrization(0, MW_SOV_SYM_LB_SS, SSHT_DL_TRAPANI),
      parametrization(0, MW_SOV_SYM_SS_POLE, SSHT_DL_RISBO),
      parametrization(0, MW_SOV_SYM_SS_POLE, SSHT_DL_TRAPANI),

      parametrization(1, DH_SOV, SSHT_DL_RISBO),
      parametrization(1, GL_SOV, SSHT_DL_RISBO),
      parametrization(1, MW_SOV_SYM, SSHT_DL_RISBO),
      parametrization(1, MW_SOV_SYM, SSHT_DL_TRAPANI),
      parametrization(1, MW_SOV_SYM_LB, SSHT_DL_RISBO),
      parametrization(1, MW_SOV_SYM_LB, SSHT_DL_TRAPANI),
      parametrization(1, MW_SOV_SYM_POLE, SSHT_DL_RISBO),
      parametrization(1, MW_SOV_SYM_POLE, SSHT_DL_TRAPANI),
      parametrization(1, MW_SOV_SYM_LB_SS, SSHT_DL_RISBO),
      parametrization(1, MW_SOV_SYM_LB_SS, SSHT_DL_TRAPANI),
      parametrization(1, MW_SOV_SYM_SS_POLE, SSHT_DL_RISBO),
      parametrization(1, MW_SOV_SYM_SS_POLE, SSHT_DL_TRAPANI),

      parametrization(-2, DH_SOV, SSHT_DL_RISBO),
      parametrization(-2, GL_SOV, SSHT_DL_RISBO),
      parametrization(-2, MW_SOV_SYM, SSHT_DL_RISBO),
      parametrization(-2, MW_SOV_SYM, SSHT_DL_TRAPANI),
      parametrization(-2, MW_SOV_SYM_LB, SSHT_DL_RISBO),
      parametrization(-2, MW_SOV_SYM_LB, SSHT_DL_TRAPANI),
      parametrization(-2, MW_SOV_SYM_POLE, SSHT_DL_RISBO),
      parametrization(-2, MW_SOV_SYM_POLE, SSHT_DL_TRAPANI),
      parametrization(-2, MW_SOV_SYM_LB_SS, SSHT_DL_RISBO),
      parametrization(-2, MW_SOV_SYM_LB_SS, SSHT_DL_TRAPANI),
      parametrization(-2, MW_SOV_SYM_SS_POLE, SSHT_DL_RISBO),
      parametrization(-2, MW_SOV_SYM_SS_POLE, SSHT_DL_TRAPANI),
      {NULL, NULL, NULL, NULL, NULL}};

  const int result = cmocka_run_group_tests(tests, NULL, NULL);

  const struct CMUnitTest *deletee = tests;
  while (deletee->test_func != NULL) {
    delete_parametrization(deletee);
    deletee = deletee + 1;
  }
  return result;
}
