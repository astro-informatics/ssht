

#include <stdlib.h>


#include "ssht_types.h"
#include "ssht_error.h"
#include "ssht_dl.h"

int ssht_dl_get_mmoffset(int L, ssht_dl_size_t dl_size);




double* ssht_dl_calloc(int L, ssht_dl_size_t dl_size) {
  
  double *dl;

  // Allocate memory.
  switch (dl_size) {

    case SSHT_DL_QUARTER:
      dl = (double*)calloc(L*L, sizeof(double));
      break;

    case SSHT_DL_HALF:
      dl = (double*)calloc(L*(2*L-1), sizeof(double));
      break;

    default:
      SSHT_ERROR_GENERIC("Invalid dl size") 

  }

  SSHT_ERROR_MEM_ALLOC_CHECK(dl)
  return dl;

}

int ssht_dl_get_mmoffset(int L, ssht_dl_size_t dl_size) {
    
  switch (dl_size) {     
    case SSHT_DL_QUARTER:
      return 0;

    case SSHT_DL_HALF:
      return L - 1;

    default:
      SSHT_ERROR_GENERIC("Invalid dl size") 
  }

}


int ssht_dl_get_mmstride(int L, ssht_dl_size_t dl_size) {
    
  switch (dl_size) {     
    case SSHT_DL_QUARTER:
      return L;

    case SSHT_DL_HALF:
      return 2*L - 1;

    default:
      SSHT_ERROR_GENERIC("Invalid dl size") 
  }

}



// dl = (double*)calloc(L*L, sizeof(double));
// dl[m*L + mm]
// m  in 0..(L-1)
// mm in 0..(L-1)

// dl = (double*)calloc(L*(2*L-1), sizeof(double));
// dl[m*(2*L-1) + mm + L-1], i.e. mm_offset = L-1
// m  in 0..(L-1)
// mm in -(L-1)..(L-1)

//signs[el+1]
void ssht_dl_halfpi_trapani_fill_eighth2righthalf_table(double *dl, int L,
							ssht_dl_size_t dl_size,
							int el, double *signs) {

  int m, mm, mmoff, mmstride;

  // Get mm offset and stride for accessing dl data.
  mmoff = ssht_dl_get_mmoffset(L, dl_size);
  mmstride = ssht_dl_get_mmstride(L, dl_size);

  // Diagonal symmetry to fill in quarter.
  for (m=0; m<=el; m++)
    for (mm=m+1; mm<=el; mm++)
      dl[m*mmstride + mm + mmoff] = 
	signs[m] * signs[mm] * dl[mm*mmstride + m + mmoff];


      //do m = 0, el
      //do mm = m+1, el
  //dl(m,mm) = signs(m) * signs(mm) * dl(mm,m);
  //     end do
  //  end do


  // Symmetry in mm to fill in half.
  for (mm=-el; mm<=-1; mm++)
    for (m=0; m<=el; m++)
      dl[m*mmstride + mm + mmoff] = 
	signs[el] * signs[m] * dl[m+mmstride - mm + mmoff];

  //do mm = -el, -1
  //do m = 0, el
  //dl(m,mm) = signs(el) * signs(m) * dl(m,-mm);
  //end do
  //end do


}



void ssht_dl_halfpi_trapani_fill_eight2quarter_table(double *dl, int L,
						     ssht_dl_size_t dl_size,
						     int el, double *signs) {

  
  int m, mm, mmoff, mmstride;

  // Get mm offset and stride for accessing dl data.
  mmoff = ssht_dl_get_mmoffset(L, dl_size);
  mmstride = ssht_dl_get_mmstride(L, dl_size);

  // Diagonal symmetry to fill in quarter.
  for (m=0; m<=el; m++)
    for (mm=m+1; mm<=el; mm++)
      dl[m*mmstride + mm + mmoff] = 
	signs[m] * signs[mm] * dl[mm*mmstride + m + mmoff];

}



void ssht_dl_halfpi_trapani_eighth_table(double *dl, int L, 
					 ssht_dl_size_t dl_size,
					 int el, double *sqrt_tbl) {

  int m, mm, mmoff, mmstride;
  double *dmm;
  double t1, t2;

  // Allocate temporary memory.
  dmm = (double*)calloc(el+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dmm)

  // Get mm offset and stride for accessing dl data.
  mmoff = ssht_dl_get_mmoffset(L, dl_size);
  mmstride = ssht_dl_get_mmstride(L, dl_size);

  // Compute Wigner plane.
  if (el == 0) {
    
    dl[0*mmstride + 0 + mmoff] = 1.0;

  }
  else {

    // Eqn (9) of T&N (2006).
    dmm[0] = - sqrt_tbl[2*el-1] / sqrt_tbl[2*el]
      * dl[(el-1)*mmstride + 0 + mmoff];
      //*(dl + (el-1)*L + 0); 

    // Eqn (10) of T&N (2006).
    for (mm=1; mm<=el; mm++) {
      dmm[mm] = sqrt_tbl[el] / SSHT_SQRT2 
        * sqrt_tbl[2*el-1] / sqrt_tbl[el+mm] / sqrt_tbl[el+mm-1]
	* dl[(el-1)*mmstride + (mm-1) + mmoff];
        //*(dl + (el-1)*L + mm-1);
    }

    // Initialise dl for next el.
    for (mm=0; mm<=el; mm++) {     
      dl[el*mmstride + mm + mmoff] = dmm[mm];
      //*(dl + el*L + mm) = dmm[mm];
    }

    // Eqn (11) of T&N (2006).
    for (mm=0; mm<=el; mm++) {     

      // m = el-1 case (t2 = 0). 
      m = el-1;
      dl[m*mmstride + mm + mmoff] = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
	* dl[(m+1)*mmstride + mm + mmoff];
      //*(dl + m*L + mm) = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
      //  * *(dl + (m+1)*L + mm);

      // Remaining m cases.
      for (m=el-2; m>=mm; m--) {
	t1 = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
	  * dl[(m+1)*mmstride + mm + mmoff];
	//  * *(dl + (m+1)*L + mm);
	t2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
	  * dl[(m+2)*mmstride + mm + mmoff];
	//* *(dl + (m+2)*L + mm);
	dl[m*mmstride + mm + mmoff] = t1 - t2;
	//*(dl + m*L + mm) = t1 - t2;
      }

    }

  }

  // Free memory.
  free(dmm);

}
