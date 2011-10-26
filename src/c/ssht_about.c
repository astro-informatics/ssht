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

int main(int argc, char *argv[]) {

  printf("%s\n", "==========================================================");
  printf("%s\n", "SSHT package to perform spin spherical harmonic transforms");
  printf("%s\n", "By Jason McEwen and Yves Wiaux");

  printf("%s\n", "See www.jasonmcewen.org for more information.");
  printf("%s\n", "See LICENSE.txt for license details.");

  printf("%s%s\n", "Version: ", SSHT_VERSION);
  printf("%s%s\n", "Build: ", SSHT_BUILD);
  printf("%s\n", "==========================================================");

  return 0;

}
