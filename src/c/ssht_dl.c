/*! 
 * \file ssht_dl.c
 * Functionality to compute Wigner plane and manage its storage.
 * 
 * \note 
 * To access the allocated memory correctly, attention must be
 * paid to its size.  This depends of the size of the dl plane
 * required, as specified by the ssht_dl_size_t type.  Memory is
 * allocated and accessed as follows for the following sizes.
 *  - SSHT_DL_QUARTER: 
 *     - allocated by dl = (double*)calloc(L*L, sizeof(double))
 *     - accessed by dl[m*L + mm], i.e. m in [0,(L-1)] and 
 *       mm in [0,(L-1)]
 *  - SSHT_DL_HALF: 
 *     - allocated by dl = (double*)calloc(L*(2*L-1), sizeof(double))
 *     - accessed by dl[m*(2*L-1) + mm + L-1], i.e. m in [0,(L-1)] and 
 *       mm in [-(L-1),(L-1)]
 *
 * \note The routine ssht_dl_calloc is provided to allocate space to
 * store a dl plane.  The routines ssht_dl_get_mmoffset and
 * ssht_dl_get_stride are provided to get appropriate offsets and
 * strides so that memory may always beaccess by dl[m*stride + mm +
 * offset].
 *
 * \author Jason McEwen
 */

#include <stdlib.h>
#include "ssht_types.h"
#include "ssht_error.h"
#include "ssht_dl.h"


/*!
 * Allocate memory to store dl plane.
 *
 * \note 
 * See file note regarding storage and access of dl memory.
 *
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate (specifying
 *            whether will need half or quarted dl plane, for example).
 * \retval dl Pointer to dl array with memory allocated.
 *
 * \author Jason McEwen
 */
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


/*!
 * Compute the offset used for accessing dl memory.
 *
 * \note 
 * See file note regarding storage and access of dl memory.
 *
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate (specifying
 *            whether will need half or quarted dl plane, for example).
 * \retval offset The offset to be used to access dl memory.
 *
 * \author Jason McEwen
 */
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


/*!
 * Compute the stride used for accessing dl memory.
 *
 * \note 
 * See file note regarding storage and access of dl memory.
 *
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate (specifying
 *            whether will need half or quarted dl plane, for example).
 * \retval stride The stride to be used to access dl memory.
 *
 * \author Jason McEwen
 */
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


/*!  
 * Calculates *eighth* (for m = 0:l and mm = 0:m) of lth plane of a
 * d-matrix for PI/2 using Trapani & Navaza's recursion method.  For
 * l>0, require the dl plane to be computed already with values for
 * l-1.  Also takes a table of precomputed square roots of integers to
 * avoid recomputing them.
 *
 * \param[in,out] dl Wigner plane.  On input this should be initialised
 * to the plane computed for el-1.  On output this will be replaced
 * with the computed plane for el.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] sqrt_tbl Precomputed array of square roots.  The table
 * element at index i should contain the value sqrt(i).  Values from 0
 * to 2*el+1 must be precomputed (i.e. sqrt_tbl should contian 2*el+2
 * elements).
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_dl_halfpi_trapani_eighth_table(double *dl, int L, 
					 ssht_dl_size_t dl_size,
					 int el, double *sqrt_tbl) {

  int m, mm, mmoff, mmstride;
  double *dmm;
  double t1, t2, s1, s2;

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

    // Eqn (10) of T&N (2006).
    for (mm=1; mm<=el; mm++) {
      dmm[mm] = sqrt_tbl[el] / SSHT_SQRT2 
        * sqrt_tbl[2*el-1] / sqrt_tbl[el+mm] / sqrt_tbl[el+mm-1]
	* dl[(el-1)*mmstride + (mm-1) + mmoff];
    }

    // Initialise dl for next el.
    for (mm=0; mm<=el; mm++) {     
      dl[el*mmstride + mm + mmoff] = dmm[mm];
    }

/*  LOGICAL BUT *NOT* MOST EFFICIENT ALGORITHM
    // Eqn (11) of T&N (2006).
    for (mm=0; mm<=el; mm++) {

      // m = el-1 case (t2 = 0).
      m = el-1;
      dl[m*mmstride + mm + mmoff] = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	* dl[(m+1)*mmstride + mm + mmoff];

      // Remaining m cases.
      for (m=el-2; m>=mm; m--) {
    	t1 = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	  * dl[(m+1)*mmstride + mm + mmoff];
    	t2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	  * dl[(m+2)*mmstride + mm + mmoff];
    	dl[m*mmstride + mm + mmoff] = t1 - t2;
      }

    }
*/

    // Eqn (11) of T&N (2006).
    // OPTIMISED FOR MEMORY ACCESS.
    m = el-1;
    s1 = sqrt_tbl[el-m] * sqrt_tbl[el+m+1];
    for (mm=0; mm<=el; mm++) {
      // m = el-1 case (t2 = 0).
      dl[m*mmstride + mm + mmoff] = 2e0 * mm / s1
    	* dl[(m+1)*mmstride + mm + mmoff];
    }

    for (m=el-2; m>=0; m--) {
      s1 = sqrt_tbl[el-m] * sqrt_tbl[el+m+1];
      s2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1];
      for (mm=0; mm<=m; mm++) {           
    	t1 = 2e0 * mm / s1
    	  * dl[(m+1)*mmstride + mm + mmoff];
    	t2 = s2
    	  * dl[(m+2)*mmstride + mm + mmoff];
    	dl[m*mmstride + mm + mmoff] = t1 - t2;
      }

    }

  }

  // Free memory.
  free(dmm);

}


/*!  
 * Calculates *quarter* (for m = 0:l and mm = 0:l) of lth plane of a
 * d-matrix for PI/2 using Trapani & Navaza's recursion method.  For
 * l>0, require the dl plane to be computed already with values for
 * l-1.  Also takes a table of precomputed square roots of integers to
 * avoid recomputing them.
 *
 * \param[in,out] dl Wigner plane.  On input this should be initialised
 * to the plane computed for el-1.  On output this will be replaced
 * with the computed plane for el.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] sqrt_tbl Precomputed array of square roots.  The table
 * element at index i should contain the value sqrt(i).  Values from 0
 * to 2*el+1 must be precomputed (i.e. sqrt_tbl should contian 2*el+2
 * elements).
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_dl_halfpi_trapani_quarter_table(double *dl, int L, 
					 ssht_dl_size_t dl_size,
					 int el, double *sqrt_tbl) {

  int m, mm, mmoff, mmstride;
  double *dmm;
  double t1, t2, s1, s2;

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

    // Eqn (10) of T&N (2006).
    for (mm=1; mm<=el; mm++) {
      dmm[mm] = sqrt_tbl[el] / SSHT_SQRT2 
        * sqrt_tbl[2*el-1] / sqrt_tbl[el+mm] / sqrt_tbl[el+mm-1]
	* dl[(el-1)*mmstride + (mm-1) + mmoff];
    }

    // Initialise dl for next el.
    for (mm=0; mm<=el; mm++) {     
      dl[el*mmstride + mm + mmoff] = dmm[mm];
    }

/*  LOGICAL BUT *NOT* MOST EFFICIENT ALGORITHM
    // Eqn (11) of T&N (2006).
    for (mm=0; mm<=el; mm++) {

      // m = el-1 case (t2 = 0).
      m = el-1;
      dl[m*mmstride + mm + mmoff] = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	* dl[(m+1)*mmstride + mm + mmoff];

      // Remaining m cases.
      for (m=el-2; m>=mm; m--) {
    	t1 = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	  * dl[(m+1)*mmstride + mm + mmoff];
    	t2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	  * dl[(m+2)*mmstride + mm + mmoff];
    	dl[m*mmstride + mm + mmoff] = t1 - t2;
      }

    }
*/

    // Eqn (11) of T&N (2006).
    // OPTIMISED FOR MEMORY ACCESS.
    m = el-1;
    s1 = sqrt_tbl[el-m] * sqrt_tbl[el+m+1];
    for (mm=0; mm<=el; mm++) {
      // m = el-1 case (t2 = 0).
      dl[m*mmstride + mm + mmoff] = 2e0 * mm / s1
    	* dl[(m+1)*mmstride + mm + mmoff];
    }

    for (m=el-2; m>=0; m--) {
      s1 = sqrt_tbl[el-m] * sqrt_tbl[el+m+1];
      s2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1];
      for (mm=0; mm<=el; mm++) {
    	t1 = 2e0 * mm / s1
    	  * dl[(m+1)*mmstride + mm + mmoff];
    	t2 = s2
    	  * dl[(m+2)*mmstride + mm + mmoff];
    	dl[m*mmstride + mm + mmoff] = t1 - t2;
      }

    }

  }

  // Free memory.
  free(dmm);

}


/*!
 * Fill in half Wigner plane for m = 0:l and mm = -l:l from the
 * eighth m = 0:l and mm = 0:m. Takes a table of precomputed signs to
 * avoid recomputing them.
 *
 * \param[in,out] dl Wigner plane.  On input should contain values
 * computed for one eighth of plane.  On output will contain values
 * computed for right half of plane, using symmetries to fill in these
 * values.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] signs Precomputed array of signs. The array element at
 * index i should contain the value (-1)^i.  Values from 0
 * to el must be precomputed (i.e. signs should contian el+1
 * elements).
 * \retval none
 *
 * \author Jason McEwen
 */
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

  // Symmetry in mm to fill in half.
  for (m=0; m<=el; m++)
    for (mm=-el; mm<=-1; mm++)
      dl[m*mmstride + mm + mmoff] = 
	signs[el] * signs[m] * dl[m*mmstride - mm + mmoff];

}


/*!
 * Fill in quarter Wigner plane for m = 0:l and mm = 0:l from the
 * eighth m = 0:l and mm = 0:m. Takes a table of precomputed signs to
 * avoid recomputing them.
 *
 * \param[in,out] dl Wigner plane.  On input should contain values
 * computed for one eighth of plane.  On output will contain values
 * computed for quarter of plane, using symmetries to fill in these
 * values.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] signs Precomputed array of signs. The array element at
 * index i should contain the value (-1)^i.  Values from 0
 * to el must be precomputed (i.e. signs should contian el+1
 * elements).
 * \retval none
 *
 * \author Jason McEwen
 */
void ssht_dl_halfpi_trapani_fill_eighth2quarter_table(double *dl, int L,
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


