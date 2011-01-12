
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define SQRT2 1.41421356237309504880168872420969807856967


void ssht_dl_halfpi_trapani_eighth_table(double *dl, int el, double *sqrt_tbl, int L);
void ssht_dl_halfpi_trapani_fill_eighth2quarter(double *dl, int el, int L);


int main(int argc, char *argv[]) {

  int i, el, m, mm, L;
  double *dl, *sqrt_tbl;
  time_t start, end;

  if (argc > 1) {
    L = atoi(argv[1]);
  }
  else {
    L = 128;
  }
  printf("Compute dl plane for L = %d...\n", L);

  dl = (double*)calloc(L*L, sizeof(double));
  sqrt_tbl = (double*)calloc(2*L-1, sizeof(double));

  // Compute square root table.
  for(i=0; i<2*L-1; i++) {
    *(sqrt_tbl + i) = sqrt((double)i);
  }

  // Compute quarter of dl plane recursively up to L-1
  // (compute eighth and fill to quarter using symmetry).
  time(&start);
  for(el=0; el<L; el++) {
    ssht_dl_halfpi_trapani_eighth_table(dl, el, sqrt_tbl, L);
    ssht_dl_halfpi_trapani_fill_eighth2quarter(dl, el, L);
  }
  time(&end);

  printf("Duration = %f seconds\n", difftime(end, start));


  // Print final dl plane.
  /* for(m=0; m<L; m++) { */
  /*   for(mm=0; mm<L; mm++) { */
  /*     printf("dl[%d][%d] = %e\n", m, mm, *(dl + m*L + mm)); */
  /*   } */
  /* } */

  // Free memory.
  free(dl);
  free(sqrt_tbl);    

}



void ssht_dl_halfpi_trapani_eighth_table(double *dl, int el, double *sqrt_tbl, int L) {

  int m, mm;
  double dmm[el+1];
  double t1, t2;

  if (el == 0) {
    
    *dl = 1.0;

  }
  else {

    // Eqn (9) of T&N (2006).
    dmm[0] = - sqrt_tbl[2*el-1] / sqrt_tbl[2*el]
      * *(dl + (el-1)*L + 0); 

    // Eqn (10) of T&N (2006).
    for (mm=1; mm<=el; mm++) {
      dmm[mm] = sqrt_tbl[el] / SQRT2 
	* sqrt_tbl[2*el-1] / sqrt_tbl[el+mm] / sqrt_tbl[el+mm-1]
	* *(dl + (el-1)*L + mm-1);
    }

    // Initialise dl for next el.
    for (mm=0; mm<=el; mm++) {     
      *(dl + el*L + mm) = dmm[mm];
    }

    // Eqn (11) of T&N (2006).
    for (mm=0; mm<=el; mm++) {     

      // m = el-1 case (t2 = 0). 
      m = el-1;
      *(dl + m*L + mm) = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
	* *(dl + (m+1)*L + mm);

      // Remaining m cases.
      for (m=el-2; m>=mm; m--) {
	t1 = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
	  * *(dl + (m+1)*L + mm);
	t2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
	  * *(dl + (m+2)*L + mm);
	*(dl + m*L + mm) = t1 - t2;
      }

    }

  }

}



void ssht_dl_halfpi_trapani_fill_eighth2quarter(double *dl, int el, int L) {

  int m, mm;
  double signs[L+1];

  for (m=0; m<L; m=m+2) {
    signs[m]   =  1.0;
    signs[m+1] = -1.0;
  }

  for (m=0; m<=el; m++) {
    for (mm=m+1; mm<=el; mm++) {
      *(dl + m*L + mm) = signs[m] * signs[mm] * *(dl + mm*L + m);
    }
  }

}
