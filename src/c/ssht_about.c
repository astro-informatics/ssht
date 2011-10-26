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
//#include "ssht_types.h"
#define SSHT_PROMPT ""


int main(int argc, char *argv[]) {

  printf("%s%s\n", SSHT_PROMPT, "==========================================================");
  printf("%s%s\n", SSHT_PROMPT, "SSHT package to perform spin spherical harmonic transforms");
  printf("%s%s\n", SSHT_PROMPT, "By Jason McEwen and Yves Wiaux");

  printf("%s%s\n", SSHT_PROMPT, "See www.jasonmcewen.org for more information.");
  printf("%s%s\n", SSHT_PROMPT, "See LICENSE.txt for license details.");

  printf("%s%s%s\n", SSHT_PROMPT, "Version: ", SSHT_VERSION);
  printf("%s%s%s\n", SSHT_PROMPT, "Build: ", SSHT_BUILD);
  printf("%s%s\n", SSHT_PROMPT, "==========================================================");

  return 0;

}
