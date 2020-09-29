#ifndef SSHT_INTERFACE_H
#define SSHT_INTERFACE_H

#include <stddef.h>

#include "ssht_dl.h"
#include "ssht_types.h"

/*! Identifies different implementation of the transforms. */
typedef enum {
  DH_SOV = 0,
  GL_SOV,
  MW_SOV_SYM,
  MW_SOV_SYM_POLE,
  MW_SOV_SYM_SS,
  MW_SOV_SYM_SS_POLE,
  MW_SOV_SYM_LB,
  MW_SOV_SYM_LB_SS,
} ssht_transforms;

/*! Parameters to the main ssht functions */
struct ssht_InterfaceParameters {
  /*! Maximum degree */
  int L;
  /*! Minimum degree */
  int L0;
  /*! 0 for less output */
  int verbosity;
  /*! dl plane recursion method */
  ssht_dl_method_t dl_method;
  /*! tranform method */
  ssht_transforms method;
};

/*! interface across all transforms in `ssht_tranforms`. */
void ssht_real_inverse(
    double *f,
    const ssht_complex_double *flm,
    const struct ssht_InterfaceParameters params);

/*! interface across all transforms in `ssht_tranforms`. */
void ssht_real_forward(
    ssht_complex_double *flm,
    const double *f,
    const struct ssht_InterfaceParameters params);

size_t image_space_size(const struct ssht_InterfaceParameters params);
size_t harmonic_space_size(const struct ssht_InterfaceParameters params);

void gen_flm_real_interface(
    complex double *flm, int seed, const struct ssht_InterfaceParameters params);
void method_name(char *result, const struct ssht_InterfaceParameters params);
#endif
