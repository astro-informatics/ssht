
#ifndef SSHT
#define SSHT


#ifdef __cplusplus
#define SSHT_COMPLEX(TYPE) std::complex<TYPE>
extern "C" {
#else
#include <complex>
#define SSHT_COMPLEX(TYPE) TYPE complex
#endif


#ifdef __cplusplus
}
#endif

SSHT_COMPLEX(double) twice(SSHT_COMPLEX(double) input) {
  return 2. * input;
}

#include "../../src/c/ssht_types.h"
#include "../../src/c/ssht_error.h"
#include "../../src/c/ssht_sampling.h"
#include "../../src/c/ssht_dl.h"
#include "../../src/c/ssht_core.h"
#include "../../src/c/ssht_adjoint.h"

#undef SSHT_COMPLEX
#endif
