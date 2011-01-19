



#include "ssht_types.h"


void ssht_dl_halfpi_trapani_eighth_table(double *dl, int nm, int nmm,
					 int el, double *sqrt_tbl) {

  int m, mm;
  double dmm[el+1];
  double t1, t2;

  //either: nmm = nm, or nmm = 2*nm - 1
  if (nmm == nm ) {
    mm_offset = 0;
  }
  else if (nmm == 2*nm - 1) {
    mm_offset = nm - 1;
  }
  else {
    printf("Error: ");

    /* __func__  */
    /* __FILE__ */
    /*   __LINE__ */
    /* __FUNCTION__ */
    /* __PRETTY_FUNCTION__ */

    exit(1);
  }


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
