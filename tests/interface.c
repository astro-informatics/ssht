#include <string.h>

#include "interface.h"
#include "ssht.h"
#include "utilities.h"

void ssht_real_inverse(
    double *f,
    const ssht_complex_double *flm,
    const struct ssht_InterfaceParameters params) {

  switch (params.method) {
  case DH_SOV:
    ssht_core_dh_inverse_sov_real(f, flm, params.L, params.verbosity);
    break;
  case GL_SOV:
    ssht_core_gl_inverse_sov_real(f, flm, params.L, params.verbosity);
    break;
  case MW_SOV_SYM:
    ssht_core_mw_inverse_sov_sym_real(
        f, flm, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_POLE:
    ssht_core_mw_inverse_sov_sym_real_pole(
        f + 1, f, flm, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_SS:
    ssht_core_mw_inverse_sov_sym_ss_real(
        f, flm, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_SS_POLE:
    ssht_core_mw_inverse_sov_sym_ss_real_pole(
        f + 2, f + 1, f, flm, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_LB:
    ssht_core_mw_lb_inverse_sov_sym_real(
        f, flm, params.L0, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_LB_SS:
    ssht_core_mw_lb_inverse_sov_sym_ss_real(
        f, flm, params.L0, params.L, params.dl_method, params.verbosity);
  }
}

void ssht_real_forward(
    ssht_complex_double *flm,
    const double *f,
    const struct ssht_InterfaceParameters params) {

  switch (params.method) {
  case DH_SOV:
    ssht_core_dh_forward_sov_real(flm, f, params.L, params.verbosity);
    break;
  case GL_SOV:
    ssht_core_gl_forward_sov_real(flm, f, params.L, params.verbosity);
    break;
  case MW_SOV_SYM:
    ssht_core_mw_forward_sov_conv_sym_real(
        flm, f, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_POLE:
    ssht_core_mw_forward_sov_conv_sym_real_pole(
        flm, f + 1, *f, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_SS:
    ssht_core_mw_forward_sov_conv_sym_ss_real(
        flm, f, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_SS_POLE:
    ssht_core_mw_forward_sov_conv_sym_ss_real_pole(
        flm, f + 2, *(f + 1), *f, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_LB:
    ssht_core_mw_lb_forward_sov_conv_sym_real(
        flm, f, params.L0, params.L, params.dl_method, params.verbosity);
    break;
  case MW_SOV_SYM_LB_SS:
    ssht_core_mw_lb_forward_sov_conv_sym_ss_real(
        flm, f, params.L0, params.L, params.dl_method, params.verbosity);
    break;
  }
}

size_t harmonic_space_size(const struct ssht_InterfaceParameters params) {
  switch (params.method) {
  case DH_SOV:
  case GL_SOV:
  case MW_SOV_SYM:
  case MW_SOV_SYM_POLE:
  case MW_SOV_SYM_SS:
  case MW_SOV_SYM_SS_POLE:
  case MW_SOV_SYM_LB:
  case MW_SOV_SYM_LB_SS:
    return params.L * params.L;
  }
}

size_t image_space_size(const struct ssht_InterfaceParameters params) {
  switch (params.method) {
  case DH_SOV:
    return 2 * params.L * (2 * params.L - 1);
  case GL_SOV:
  case MW_SOV_SYM:
  case MW_SOV_SYM_LB:
    return params.L * (2 * params.L - 1);
  case MW_SOV_SYM_POLE:
    return (params.L - 1) * (2 * params.L - 1) + 1;
  case MW_SOV_SYM_SS:
  case MW_SOV_SYM_LB_SS:
    return 2 * params.L * (params.L + 1);
  case MW_SOV_SYM_SS_POLE:
    return 2 * params.L * (params.L - 1) + 2;
  }
}

void gen_flm_real_interface(
    complex double *flm, int seed, const struct ssht_InterfaceParameters params) {
  switch (params.method) {
  case DH_SOV:
  case GL_SOV:
  case MW_SOV_SYM_POLE:
  case MW_SOV_SYM:
  case MW_SOV_SYM_SS:
  case MW_SOV_SYM_SS_POLE:
    gen_flm_real(flm, params.L, seed);
    break;
  case MW_SOV_SYM_LB:
  case MW_SOV_SYM_LB_SS:
    gen_lb_flm_real(flm, params.L0, params.L, seed);
    break;
  }
}

void method_name(char *result, const struct ssht_InterfaceParameters params) {
  char dl[30];
  switch (params.dl_method) {
  case SSHT_DL_RISBO:
    strcpy(dl, "(with Risbo d-l plane)");
    break;
  case SSHT_DL_TRAPANI:
    strcpy(dl, "(with Trapani d-l plane)");
    break;
  }
  switch (params.method) {
  case DH_SOV:
    strcpy(result, "DH");
    break;
  case GL_SOV:
    strcpy(result, "GH");
    break;
  case MW_SOV_SYM_POLE:
    sprintf(result, "Symmetric MW with poles%s", dl);
    break;
  case MW_SOV_SYM:
    sprintf(result, "Symmetric MW%s", dl);
    break;
  case MW_SOV_SYM_SS_POLE:
    sprintf(result, "MW including all symmetries and poles%s", dl);
    break;
  case MW_SOV_SYM_LB:
    sprintf(result, "Symmetric MW with lower bound%s", dl);
    break;
  case MW_SOV_SYM_SS:
    sprintf(result, "MW including all symmetries%s", dl);
    break;
  case MW_SOV_SYM_LB_SS:
    sprintf(result, "MW including all symmetries and lower bound%s", dl);
    break;
  }
}
