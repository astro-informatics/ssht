// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details

/*! 
 * \file ssht_about.c
 * Print information about the SSHT package, including version
 * and build numbers. 
 *
 * Usage: ssht_about
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#include <stdio.h>
#ifdef BUILT_WITH_CMAKE
#include "ssht/ssht_version.h"
#endif

int main(void) {

#ifndef BUILT_WITH_CMAKE
  printf("%s\n", "==========================================================");
  printf("%s\n", "SSHT package to perform spin spherical harmonic transforms");
  printf("%s\n", "By Jason McEwen, Chris Wallis, Martin Buttner, ");
  printf("%s\n", "Boris Leistedt, and Yves Wiaux");

  printf("%s\n", "See www.jasonmcewen.org for more information.");
  printf("%s\n", "See LICENSE.txt for license details.");

  printf("%s%s\n", "Version: ", SSHT_VERSION);
  printf("%s%s\n", "Build: ", SSHT_BUILD);
  printf("%s\n", "==========================================================");
#else
  printf("%s", ssht_info());
#endif

  return 0;

}
