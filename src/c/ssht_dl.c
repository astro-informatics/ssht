


#include "ssht_types.h"
#include "ssht_dl.h"


// dl = (double*)calloc(L*L, sizeof(double));
// dl[m*L + mm]
// m  in 0..(L-1)
// mm in 0..(L-1)

// dl = (double*)calloc(L*(2*L-1), sizeof(double));
// dl[m*(2*L-1) + mm + L-1], i.e. mm_offset = L-1
// m  in 0..(L-1)
// mm in -(L-1)..(L-1)




double* ssht_dl_calloc(int L, ssht_dl_size_t dl_size) {
  
  double *dl = NULL;

  // Allocate memory.
  switch (dl_size) {

    case SSHT_DL_QUARTER:
      dl = (double*)calloc(L*L, sizeof(double));
      break;

    case SSHT_DL_HALF:
      dl = (double*)calloc(L*(2*L-1), sizeof(double));
      break;

    default:
      printf("Error:");
      exit(1);

  }

  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  return dl;

}



void ssht_dl_halfpi_trapani_eighth_table(double *dl, int L, 
					 ssht_dl_size_t dl_size,
					 int el, double *sqrt_tbl) {

  int m, mm;
  double *dmm;
  double t1, t2;

  // Allocate temporary memory.
  dmm = (double*)calloc(el+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dmm)



  //either: nmm = nm, or nmm = 2*nm - 1
  if (nmm == nm) {
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
    
    dl[0*] = 1.0;

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

  // Free memory.
  free(dmm);

}
