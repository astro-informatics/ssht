// SSHT package to perform spin spherical harmonic transforms
// Copyright (C) 2011  Jason McEwen
// See LICENSE.txt for license details


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
 *  - SSHT_DL_QUARTER_EXTENDED: 
 *     - allocated by dl = (double*)calloc((L+1)*(L+3), sizeof(double))
 *     - require extended size when using Risbo recursion for
 *       partial Wigner plane
 *     - accessed by dl[(m + (L-1))*(L+3) + mm + L-1], i.e. m in 
 *       [-(L-1),0] and mm in [-(L-1),0]
 *  - SSHT_DL_HALF: 
 *     - allocated by dl = (double*)calloc(L*(2*L-1), sizeof(double))
 *     - accessed by dl[m*(2*L-1) + mm + L-1], i.e. m in [0,(L-1)] and 
 *       mm in [-(L-1),(L-1)]
 *  - SSHT_DL_FULL: 
 *     - allocated by dl = (double*)calloc((2*L-1)*(2*L-1), sizeof(double))
 *     - accessed by dl[(m + (L-1))*(2*L-1) + mm + L-1], i.e. 
 *       m in [-(L-1),(L-1)] and mm in [-(L-1),(L-1)]
 *
 * \note The routine ssht_dl_calloc is provided to allocate space to
 * store a dl plane.  The routines ssht_dl_get_offset and
 * ssht_dl_get_stride are provided to get appropriate offsets and
 * strides.
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */

#include <stdlib.h>
#include <math.h>
#include "ssht/ssht_types.h"
#include "ssht/ssht_error.h"
#include "ssht/ssht_dl.h"

double logfact(int n);


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
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
double* ssht_dl_calloc(int L, ssht_dl_size_t dl_size) {
  
  double *dl;

  // Allocate memory.
  switch (dl_size) {

    case SSHT_DL_QUARTER:
      dl = (double*)calloc(L*L, sizeof(double));
      break;

    case SSHT_DL_QUARTER_EXTENDED:
      dl = (double*)calloc((L+2)*(L+2), sizeof(double));
      break;

    case SSHT_DL_HALF:
      dl = (double*)calloc(L*(2*L-1), sizeof(double));
      break;

    case SSHT_DL_FULL:
      dl = (double*)calloc((2*L-1)*(2*L-1), sizeof(double));
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
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_dl_get_offset(int L, ssht_dl_size_t dl_size) {
    
  switch (dl_size) {     
    case SSHT_DL_QUARTER:
      return 0;

    case SSHT_DL_QUARTER_EXTENDED:
      return L - 1;

    case SSHT_DL_HALF:
      return L - 1;

    case SSHT_DL_FULL:
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
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
int ssht_dl_get_stride(int L, ssht_dl_size_t dl_size) {
    
  switch (dl_size) {     
    case SSHT_DL_QUARTER:
      return L;

    case SSHT_DL_QUARTER_EXTENDED:
      return L + 2;

    case SSHT_DL_HALF:
      return 2*L - 1;

    case SSHT_DL_FULL:
      return 2*L - 1;

    default:
      SSHT_ERROR_GENERIC("Invalid dl size") 
  }

}


/*!  
 * Calculates (for m = -l:l and mm = -l:l) lth plane of a
 * d-matrix for argument beta using Risbo's recursion method.  For
 * l>0, require the dl plane to be computed already with values for
 * l-1.  Also takes a table of precomputed square roots of integers to
 * avoid recomputing them.
 *
 * \param[in,out] dl Wigner plane.  On input this should be initialised
 * to the plane computed for el-1.  On output this will be replaced
 * with the computed plane for el.
 * \param[in] beta Angle to compute Wigner line for.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] sqrt_tbl Precomputed array of square roots.  The table
 * element at index i should contain the value sqrt(i).  Values from 0
 * to 2*el must be precomputed (i.e. sqrt_tbl should contian 2*el+1
 * elements).
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_beta_risbo_full_table(double *dl, double beta, int L, 
				      ssht_dl_size_t dl_size,
				      int el, double *sqrt_tbl) {

  int offset, stride;
  double cosb, sinb, coshb, sinhb;
  int i, j, k;
  double rj, dlj, ddj;
  double *dd;
    
  // Get mm offset and stride for accessing dl data.
  offset = ssht_dl_get_offset(L, dl_size);
  stride = ssht_dl_get_stride(L, dl_size);

  // Compute Wigner plane.
  if (el == 0) {
    
    dl[(0+offset)*stride + 0 + offset] = 1.0;

  }
  else if (el == 1) {

    cosb = cos(beta);
    sinb = sin(beta);
    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);
   
    dl[(-1 + offset)*stride - 1 + offset] = coshb * coshb;
    dl[(-1 + offset)*stride + 0 + offset] = sinb / SSHT_SQRT2;
    dl[(-1 + offset)*stride + 1 + offset] = sinhb * sinhb;

    dl[(0 + offset)*stride - 1 + offset] = -sinb / SSHT_SQRT2;
    dl[(0 + offset)*stride + 0 + offset] = cosb;
    dl[(0 + offset)*stride + 1 + offset] = sinb / SSHT_SQRT2;

    dl[(1 + offset)*stride - 1 + offset] = sinhb * sinhb;
    dl[(1 + offset)*stride + 0 + offset] = -sinb / SSHT_SQRT2;
    dl[(1 + offset)*stride + 1 + offset] = coshb * coshb;

  }
  else {

    coshb = -cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    // Initialise the plane of the dl-matrix to 0.0 for the recursion
    // from l - 1 to l - 1/2.
    dd = (double*)calloc((2*el+2)*(2*el+2), sizeof(double));
    SSHT_ERROR_MEM_ALLOC_CHECK(dd)

    j = 2*el - 1;
    rj = (double) j;
    for (k=0; k<=j-1; k++) {
      for (i=0; i<=j-1; i++) {
	dlj = dl[(k-(el-1)+offset)*stride + i-(el-1) + offset] / rj;
	dd[k*(2*el+2) + i] +=
	  sqrt_tbl[j-i] * sqrt_tbl[j-k] * dlj * coshb;
	dd[k*(2*el+2) + i+1] -=
	  sqrt_tbl[i+1] * sqrt_tbl[j-k] * dlj * sinhb;
	dd[(k+1)*(2*el+2) + i] +=
	  sqrt_tbl[j-i] * sqrt_tbl[k+1] * dlj * sinhb;
	dd[(k+1)*(2*el+2) + i+1] +=
	  sqrt_tbl[i+1] * sqrt_tbl[k+1] * dlj * coshb;
      }
    }

    // Having constructed the d^(l+1/2) matrix in dd, do the second
    // half-step recursion from dd to dl. Start by initilalising  
    // the plane of the dl-matrix to 0.0.
    for (k=-el; k<=el; k++) 
      for (i=-el; i<=el; i++)
	dl[(k+offset)*stride + i + offset] = 0.0;

    j = 2*el;
    rj = (double) j;
    for (k=0; k<=j-1; k++) {
      for (i=0; i<=j-1; i++) {
	ddj = dd[k*(2*el+2) + i] / rj;
	dl[(k-el+offset)*stride + i-el + offset] +=
	  sqrt_tbl[j-i] * sqrt_tbl[j-k] * ddj * coshb;
	dl[(k-el+offset)*stride + i+1-el + offset] -=
	  sqrt_tbl[i+1] * sqrt_tbl[j-k] * ddj * sinhb;
	dl[(k+1-el+offset)*stride + i-el + offset] +=
	  sqrt_tbl[j-i] * sqrt_tbl[k+1] * ddj * sinhb;
	dl[(k+1-el+offset)*stride + i+1-el + offset] +=
	  sqrt_tbl[i+1] * sqrt_tbl[k+1] * ddj * coshb;
      }
    }

    // Free temporary memory.
    free(dd);
  }

}


/*!  
 * Calculates (for m = -l:l and mm = -l:l) lth plane of a d-matrix for
 * argument beta using Risbo's recursion method.  Only half of the
 * plane is computed by recusion and symmetry is used to fill the
 * remaining plane.  For l>0, require the dl plane to be computed
 * already with values for l-1.  Also takes a table of precomputed
 * square roots of integers and signs to avoid recomputing them.
 *
 * \param[in,out] dl Wigner plane.  On input this should be initialised
 * to the plane computed for el-1.  On output this will be replaced
 * with the computed plane for el.
 * \param[in] beta Angle to compute Wigner line for.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] sqrt_tbl Precomputed array of square roots.  The table
 * element at index i should contain the value sqrt(i).  Values from 0
 * to 2*el must be precomputed (i.e. sqrt_tbl should contian 2*el+1
 * elements).
 * \param[in] signs Precomputed array of signs. The array element at
 * index i should contain the value (-1)^i.  Values from 0
 * to el must be precomputed (i.e. signs should contian el+1
 * elements).
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_beta_risbo_half_table(double *dl, double beta, int L, 
				   ssht_dl_size_t dl_size,
				   int el, double *sqrt_tbl,
				   double *signs) {

  int offset, stride;
  double cosb, sinb, coshb, sinhb;
  int i, j, k;
  double rj, dlj, ddj;
  double *dd;
  int m, mm;

  // Get mm offset and stride for accessing dl data.
  offset = ssht_dl_get_offset(L, dl_size);
  stride = ssht_dl_get_stride(L, dl_size);

  // Compute Wigner plane.
  if (el == 0) {
    
    dl[(0+offset)*stride + 0 + offset] = 1.0;

  }
  else if (el == 1) {

    cosb = cos(beta);
    sinb = sin(beta);
    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);
   
    dl[(-1 + offset)*stride - 1 + offset] = coshb * coshb;
    dl[(-1 + offset)*stride + 0 + offset] = sinb / SSHT_SQRT2;
    dl[(-1 + offset)*stride + 1 + offset] = sinhb * sinhb;

    dl[(0 + offset)*stride - 1 + offset] = -sinb / SSHT_SQRT2;
    dl[(0 + offset)*stride + 0 + offset] = cosb;
    dl[(0 + offset)*stride + 1 + offset] = sinb / SSHT_SQRT2;

    dl[(1 + offset)*stride - 1 + offset] = sinhb * sinhb;
    dl[(1 + offset)*stride + 0 + offset] = -sinb / SSHT_SQRT2;
    dl[(1 + offset)*stride + 1 + offset] = coshb * coshb;

  }
  else {

    coshb = -cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    // Initialise the plane of the dl-matrix to 0.0 for the recursion
    // from l - 1 to l - 1/2.
    dd = (double*)calloc((2*el+2)*(2*el+2), sizeof(double));
    SSHT_ERROR_MEM_ALLOC_CHECK(dd)

    j = 2*el - 1;
    rj = (double) j;
    for (k=0; k<=j-1; k++) {
      for (i=0; i<=el; i++) {
	dlj = dl[(k-(el-1)+offset)*stride + i-(el-1) + offset] / rj;
	dd[k*(2*el+2) + i] +=
	  sqrt_tbl[j-i] * sqrt_tbl[j-k] * dlj * coshb;
	dd[k*(2*el+2) + i+1] -=
	  sqrt_tbl[i+1] * sqrt_tbl[j-k] * dlj * sinhb;
	dd[(k+1)*(2*el+2) + i] +=
	  sqrt_tbl[j-i] * sqrt_tbl[k+1] * dlj * sinhb;
	dd[(k+1)*(2*el+2) + i+1] +=
	  sqrt_tbl[i+1] * sqrt_tbl[k+1] * dlj * coshb;
      }
    }

    // Having constructed the d^(l+1/2) matrix in dd, do the second
    // half-step recursion from dd to dl. Start by initilalising  
    // the plane of the dl-matrix to 0.0.
    for (k=-el; k<=el; k++) 
      for (i=-el; i<=el; i++)
	dl[(k+offset)*stride + i + offset] = 0.0;

    j = 2*el;
    rj = (double) j;
    for (k=0; k<=j-1; k++) {
      for (i=0; i<=el; i++) {
	ddj = dd[k*(2*el+2) + i] / rj;
	dl[(k-el+offset)*stride + i-el + offset] +=
	  sqrt_tbl[j-i] * sqrt_tbl[j-k] * ddj * coshb;
	dl[(k-el+offset)*stride + i+1-el + offset] -=
	  sqrt_tbl[i+1] * sqrt_tbl[j-k] * ddj * sinhb;
	dl[(k+1-el+offset)*stride + i-el + offset] +=
	  sqrt_tbl[j-i] * sqrt_tbl[k+1] * ddj * sinhb;
	dl[(k+1-el+offset)*stride + i+1-el + offset] +=
	  sqrt_tbl[i+1] * sqrt_tbl[k+1] * ddj * coshb;
      }
    }

    // Fill top half of plane using symmetry.
    for (m=-el; m<=el; m++) 
      for (mm=1; mm<=el; mm++) 
        dl[(m+offset)*stride + mm + offset] = 
          signs[abs(m)] * signs[abs(mm)]
          * dl[(-m+offset)*stride - mm + offset];

    // Free temporary memory.
    free(dd);
  }

}


/*!  
 * Calculates (for m = -l:0 and mm = -l:m) lth plane of a d-matrix for
 * argument beta using Risbo's recursion method.  For l>0, require the
 * dl plane to be computed already with values for l-1.  Also takes a
 * table of precomputed square roots of integers and signs to avoid
 * recomputing them.
 *
 * \note Risbo's recursion requires the boundaries of the eighth of
 * the Wigner plane to be extended slightly.  Consequently, this
 * routine should only be used with dl_size either
 * SSHT_DL_QUARTER_EXTENDED or SSHT_DL_FULL.
 *
 * \param[in,out] dl Wigner plane.  On input this should be initialised
 * to the plane computed for el-1.  On output this will be replaced
 * with the computed plane for el.
 * \param[in] beta Angle to compute Wigner line for.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] sqrt_tbl Precomputed array of square roots.  The table
 * element at index i should contain the value sqrt(i).  Values from 0
 * to 2*el must be precomputed (i.e. sqrt_tbl should contian 2*el+1
 * elements).
 * \param[in] signs Precomputed array of signs. The array element at
 * index i should contain the value (-1)^i.  Values from 0
 * to el must be precomputed (i.e. signs should contian el+1
 * elements).
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_beta_risbo_eighth_table(double *dl, double beta, int L, 
				      ssht_dl_size_t dl_size,
				      int el, double *sqrt_tbl,
				      double *signs) {

  int offset, stride, imax;
  double cosb, sinb, coshb, sinhb;
  int i, j, k;
  double rj, dlj, ddj;
  double *dd;
  int m, mm;

  // Get mm offset and stride for accessing dl data.
  offset = ssht_dl_get_offset(L, dl_size);
  stride = ssht_dl_get_stride(L, dl_size);

  // Compute Wigner plane.
  if (el == 0) {
    
    dl[(0+offset)*stride + 0 + offset] = 1.0;

  }
  else if (el == 1) {

    cosb = cos(beta);
    sinb = sin(beta);
    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);
   
    dl[(-1 + offset)*stride - 1 + offset] = coshb * coshb;
    dl[(-1 + offset)*stride + 0 + offset] = sinb / SSHT_SQRT2;
    dl[(-1 + offset)*stride + 1 + offset] = sinhb * sinhb;

    dl[(0 + offset)*stride - 1 + offset] = -sinb / SSHT_SQRT2;
    dl[(0 + offset)*stride + 0 + offset] = cosb;
    dl[(0 + offset)*stride + 1 + offset] = sinb / SSHT_SQRT2;

    dl[(1 + offset)*stride - 1 + offset] = sinhb * sinhb;
    dl[(1 + offset)*stride + 0 + offset] = -sinb / SSHT_SQRT2;
    dl[(1 + offset)*stride + 1 + offset] = coshb * coshb;

  }
  else {

    coshb = -cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    // Initialise the plane of the dl-matrix to 0.0 for the recursion
    // from l - 1 to l - 1/2.
    dd = (double*)calloc((el+3)*(el+3), sizeof(double));
    SSHT_ERROR_MEM_ALLOC_CHECK(dd)

    j = 2*el - 1;
    rj = (double) j;
    for (k=0; k<=el; k++) {
      if (k==el)
	imax = k+1;
      else
	imax = k+2;
      for (i=0; i<=imax; i++) {
	dlj = dl[(k-(el-1)+offset)*stride + i-(el-1) + offset] / rj;	
	dd[k*(el+3) + i] +=
	  sqrt_tbl[j-i] * sqrt_tbl[j-k] * dlj * coshb;
	dd[k*(el+3) + i+1] -=
	  sqrt_tbl[i+1] * sqrt_tbl[j-k] * dlj * sinhb;
	dd[(k+1)*(el+3) + i] +=
	  sqrt_tbl[j-i] * sqrt_tbl[k+1] * dlj * sinhb;
	dd[(k+1)*(el+3) + i+1] +=
	  sqrt_tbl[i+1] * sqrt_tbl[k+1] * dlj * coshb;
      }
    }

    // Having constructed the d^(l+1/2) matrix in dd, do the second
    // half-step recursion from dd to dl. Start by initilalising  
    // the plane of the dl-matrix to 0.0.
    for (k=-el; k<=1; k++) 
      for (i=-el; i<=3; i++)
	dl[(k+offset)*stride + i + offset] = 0.0;

    j = 2*el;
    rj = (double) j;
    for (k=0; k<=el; k++) {
      for (i=0; i<=k+1; i++) {
	ddj = dd[k*(el+3) + i] / rj;
	dl[(k-el+offset)*stride + i-el + offset] +=
	  sqrt_tbl[j-i] * sqrt_tbl[j-k] * ddj * coshb;
	dl[(k-el+offset)*stride + i+1-el + offset] -=
	  sqrt_tbl[i+1] * sqrt_tbl[j-k] * ddj * sinhb;
	dl[(k+1-el+offset)*stride + i-el + offset] +=
	  sqrt_tbl[j-i] * sqrt_tbl[k+1] * ddj * sinhb;
	dl[(k+1-el+offset)*stride + i+1-el + offset] +=
	  sqrt_tbl[i+1] * sqrt_tbl[k+1] * ddj * coshb;
      }
    }

    // Extend dl plane about boundaries, using symmetries, to the
    // extents required for recursion.

    // Extend above diagonal by two in mm.
    for (m=-el; m<=0; m++)
      for (mm=m+1; mm<=m+2; mm++)
    	dl[(m+offset)*stride + mm + offset] =
	  signs[abs(m)] * signs[abs(mm)] 
	  * dl[(mm+offset)*stride + m + offset];

    // Extend right by one in m.
    for (m=1; m<=1; m++)
      for (mm=-el; mm<=0; mm++)
    	dl[(m+offset)*stride + mm + offset] =
	  signs[abs(el)] * signs[abs(mm)] 
	  * dl[(-m+offset)*stride + mm + offset];

    // Extend up by one in mm.
    for (m=-el; m<=1; m++)
      for (mm=1; mm<=1; mm++)
    	dl[(m+offset)*stride + mm + offset] =
    	  signs[abs(el)] * signs[abs(m)]
	  * dl[(m+offset)*stride - mm + offset];

    // Free temporary memory.
    free(dd);
  }

}


/*!
 * Fill in quarter Wigner plane for m = 0:l and mm = 0:l from the
 * eighth m = -l:0 and mm = -l:m. Takes a table of precomputed signs to
 * avoid recomputing them.
 *
 * \note It is necessary to fill the quarter Wigner plane in separate
 * memory since the eighth plane for m = -l:0 and mm = -l:m with be
 * required in subsequent recursions.
 *
 * \param[out] dl Quarter Wigner plane, computed by copying dl8 and
 * extending with symmetries.
 * \param[in] dl8 Eighth Wigner plane.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the dl memory.
 * \param[in] dl8_size Size type of the dl8 memory.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] signs Precomputed array of signs. The array element at
 * index i should contain the value (-1)^i.  Values from 0
 * to el must be precomputed (i.e. signs should contian el+1
 * elements).
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_beta_risbo_fill_eighth2quarter_table(double *dl, 
						  double *dl8,
						  int L,
						  ssht_dl_size_t dl_size,
						  ssht_dl_size_t dl8_size,
						  int el, 
						  double *signs) {

  int offset, stride, offset8, stride8;
  int m, mm;

  // Get mm offset and stride for accessing dl data.
  offset = ssht_dl_get_offset(L, dl_size);
  stride = ssht_dl_get_stride(L, dl_size);
  offset8 = ssht_dl_get_offset(L, dl8_size);
  stride8 = ssht_dl_get_stride(L, dl8_size);

  // Symmetry through origin to get eighth of the required quarter
  // (first quadrant).
  for (m=0; m<=el; m++)
    for (mm=m; mm<=el; mm++)
      dl[(m+offset)*stride + mm + offset] =
	signs[m] * signs[mm]
  	* dl8[(-m+offset8)*stride8 - mm + offset8];

  // Diagonal symmetry to fill remaining quarter.
  for (m=0; m<=el; m++)
    for (mm=0; mm<=m-1; mm++)
      dl[(m+offset)*stride + mm + offset] =
	signs[m] * signs[mm]
  	* dl[(mm+offset)*stride + m + offset];

}


/*!  
 * Calculates (for m = -l:l and mm = -l:l) lth plane of a d-matrix for
 * argument beta using the recursion method given in Kostelec and
 * Rockmore (2010) (see equations (4.5)-(4.9)).  For l>1, require the
 * dl plane to be computed already with values for l-1 *and* l-2.
 * Also takes a table of precomputed square roots of integers and
 * signs to avoid recomputing them.
 *
 * \param[in,out] dlm1p1 Wigner plane.  On input this should be initialised
 * to the plane computed for el-2.  On output this will be replaced
 * with the computed plane for el.
 * \param[in] dl Wigner plane already computed for el-1.
 * \param[in] beta Angle to compute Wigner plan for.
 * \param[in] L Harmonic band-limit.
 * \param[in] dl_size Size type of the memory to allocate.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] sqrt_tbl Precomputed array of square roots.  The table
 * element at index i should contain the value sqrt(i).  Values from 0
 * to 2*el must be precomputed (i.e. sqrt_tbl should contian 2*el+1
 * elements).
 * \param[in] signs Precomputed array of signs. The array element at
 * index i should contain the value (-1)^i.  Values from 0
 * to el must be precomputed (i.e. signs should contian el+1
 * elements).
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_beta_kostelec_full_table(double *dlm1p1, double *dl, 
				      double beta, int L, 
				      ssht_dl_size_t dl_size,
				      int el, 
				      double *sqrt_tbl, double *signs) {

  int offset, stride;
  double cosb, sinb, coshb, sinhb;
  double lnAlm, lnfact2el;
  int m, mm, elm1;
  double elr, elm1r;
    
  // Get mm offset and stride for accessing dl data.
  offset = ssht_dl_get_offset(L, dl_size);
  stride = ssht_dl_get_stride(L, dl_size);

  // Compute Wigner plane.
  if (el == 0) {
    
    dlm1p1[(0+offset)*stride + 0 + offset] = 1.0;

  }
  else if (el == 1) {

    cosb = cos(beta);
    sinb = sin(beta);
    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    // These terms follow directly from the boundary conditions and
    // one recursion to get the central term.
   
    dlm1p1[(-1 + offset)*stride - 1 + offset] = coshb * coshb;
    dlm1p1[(-1 + offset)*stride + 0 + offset] = sinb / SSHT_SQRT2;
    dlm1p1[(-1 + offset)*stride + 1 + offset] = sinhb * sinhb;

    dlm1p1[(0 + offset)*stride - 1 + offset] = -sinb / SSHT_SQRT2;
    dlm1p1[(0 + offset)*stride + 0 + offset] = cosb;
    dlm1p1[(0 + offset)*stride + 1 + offset] = sinb / SSHT_SQRT2;

    dlm1p1[(1 + offset)*stride - 1 + offset] = sinhb * sinhb;
    dlm1p1[(1 + offset)*stride + 0 + offset] = -sinb / SSHT_SQRT2;
    dlm1p1[(1 + offset)*stride + 1 + offset] = coshb * coshb;

  }
  else {

    cosb = cos(beta);
    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    elr = (double) el;
    elm1 = el - 1;
    elm1r = (double) elm1;

    // Recurse over all of plane except boundaries.
    for (m=-(el-1); m<=el-1; m++) {
      for (mm=-(el-1); mm<=el-1; mm++) {

	// Compute 3-term recursion.
	dlm1p1[(m + offset)*stride + mm + offset] = 
	  (cosb - m*mm/(elm1r*elr)) * dl[(m + offset)*stride + mm + offset] 
	  - 
	  sqrt_tbl[elm1+m] * sqrt_tbl[elm1-m] * sqrt_tbl[elm1+mm] * sqrt_tbl[elm1-mm] 
	  / (elm1r * (2.0*elm1r + 1.0))
	  * dlm1p1[(m + offset)*stride + mm + offset];

	// Perform scaling.
	dlm1p1[(m + offset)*stride + mm + offset] *= 
	  el * (2*elm1 + 1.0)
	  / (sqrt_tbl[el-m] * sqrt_tbl[el+m] * sqrt_tbl[el-mm] * sqrt_tbl[el+mm]);

      }
    }

    // Compute boundaries.
    lnfact2el = logfact(2*el);
    for (m=-el; m<=el; m++) {

      lnAlm = (lnfact2el - logfact(el+m) - logfact(el-m)) / 2.0;

      // Right line.
      dlm1p1[(el + offset)*stride + m + offset] =  
	signs[el] * signs[abs(m)]
	* exp(lnAlm + (el+m)*log(coshb) + (el-m)*log(sinhb));

      // Left line.
      dlm1p1[(-el + offset)*stride + m + offset] =  
	exp(lnAlm + (el-m)*log(coshb) + (el+m)*log(sinhb));

      // Top line.
      dlm1p1[(m + offset)*stride + el + offset] =  
	exp(lnAlm + (el+m)*log(coshb) + (el-m)*log(sinhb));

      // Bottom line.
      dlm1p1[(m + offset)*stride - el + offset] =  
	signs[el] * signs[abs(m)]
	* exp(lnAlm + (el-m)*log(coshb) + (el+m)*log(sinhb));
    }

  }

}


/*!  
 * Calculates line of d-matrix for all m = -l:l and given mm = -l:l for 
 * argument beta using the recursion method given in Kostelec and
 * Rockmore (2010) (see equations (4.5)-(4.9)).  For l>abs(mm), require the
 * dl plane to be computed already with values for l-1 *and* l-2.
 * Also takes a table of precomputed square roots of integers and
 * signs to avoid recomputing them.
 *
 * \param[in,out] dlm1p1_line Wigner line.  On input this should be initialised
 * to the line computed for el-2.  On output this will be replaced
 * with the computed line for el.
 * \param[in] dl Wigner plane already computed for el-1.  
 * \param[in] beta Angle to compute Wigner line for.
 * \param[in] L Harmonic band-limit.
 * \param[in] mm Azimuthal harmonic index to compute Wigner line for.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] sqrt_tbl Precomputed array of square roots.  The table
 * element at index i should contain the value sqrt(i).  Values from 0
 * to 2*el must be precomputed (i.e. sqrt_tbl should contian 2*el+1
 * elements).
 * \param[in] signs Precomputed array of signs. The array element at
 * index i should contain the value (-1)^i.  Values from 0
 * to el must be precomputed (i.e. signs should contian el+1
 * elements).
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_beta_kostelec_line_table(double *dlm1p1_line, double *dl_line, 
				      double beta, int L, int mm, int el, 
				      double *sqrt_tbl, double *signs) {

  int offset;
  double cosb, sinb, coshb, sinhb;
  double lnAlm, lnAlmm, lnfact2el;
  int m, elm1;
  double elr, elm1r;
  
  // Compute m offset for accessing dl line.
  offset = L-1;

  // Compute Wigner plane.
  if (el < abs(mm)) {
    // Do nothing (dl line should remain zero).
    return;
  }
  else if (el == 1) {

    cosb = cos(beta);
    sinb = sin(beta);
    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    if (mm == -1) {
      dlm1p1_line[-1 + offset] = coshb * coshb;
      dlm1p1_line[ 0 + offset] = -sinb / SSHT_SQRT2;
      dlm1p1_line[ 1 + offset] = sinhb * sinhb;
    }
    else if (mm == 0) {
      dlm1p1_line[-1 + offset] = sinb / SSHT_SQRT2;
      dlm1p1_line[ 0 + offset] = cosb;
      dlm1p1_line[ 1 + offset] = -sinb / SSHT_SQRT2;
    }
    else {
      // mm == +1
      dlm1p1_line[-1 + offset] = sinhb * sinhb;
      dlm1p1_line[ 0 + offset] = sinb / SSHT_SQRT2;
      dlm1p1_line[ 1 + offset] = coshb * coshb;
    }

  }
  else if (el == abs(mm)) {

    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    // Initalise line using equation (4.8) or (4.9) from K&R (2010).
    if (mm >= 0) {

      // Initialise using equation (4.8), i.e. top line.
      lnfact2el = logfact(2*el);
      for (m=-el; m<=el; m++) {
	lnAlm = (lnfact2el - logfact(el+m) - logfact(el-m)) / 2.0;
	dlm1p1_line[m + offset] = 
	  exp(lnAlm + (el+m)*log(coshb) + (el-m)*log(sinhb));
      }

    }
    else {

      // Initialise using equation (4.9), i.e. bottom line.
      lnfact2el = logfact(2*el);
      for (m=-el; m<=el; m++) {
	lnAlm = (lnfact2el - logfact(el+m) - logfact(el-m)) / 2.0;
	dlm1p1_line[m + offset] = 
	  signs[el] * signs[abs(m)]
	  * exp(lnAlm + (el-m)*log(coshb) + (el+m)*log(sinhb));
      }

    }

  }
  else {

    cosb = cos(beta);
    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    elr = (double) el;
    elm1 = el - 1;
    elm1r = (double) elm1;

    // Recuse over line.
    for (m=-(el-1); m<=el-1; m++) {

      // Compute 3-term recursion.
      dlm1p1_line[m + offset] = 
	(cosb - m*mm/(elm1r*elr)) * dl_line[m + offset] 
	- 
	sqrt_tbl[elm1+m] * sqrt_tbl[elm1-m] * sqrt_tbl[elm1+mm] * sqrt_tbl[elm1-mm] 
	/ (elm1r * (2.0*elm1r + 1.0))
	* dlm1p1_line[m + offset];
      
      // Perform scaling.
      dlm1p1_line[m + offset] *= 
	el * (2*elm1 + 1.0)
	/ (sqrt_tbl[el-m] * sqrt_tbl[el+m] * sqrt_tbl[el-mm] * sqrt_tbl[el+mm]);

    }

    // Compute edges...
    lnfact2el = logfact(2*el);
    lnAlmm = (lnfact2el - logfact(el+mm) - logfact(el-mm)) / 2.0;
    
    // Left edge.
    dlm1p1_line[-el + offset] =  
      exp(lnAlmm + (el-mm)*log(coshb) + (el+mm)*log(sinhb));

    // Right edge.
    dlm1p1_line[el + offset] =  
      signs[el] * signs[abs(mm)]
      * exp(lnAlmm + (el+mm)*log(coshb) + (el-mm)*log(sinhb));

  }

}


/*!  
 * Calculates half line of d-matrix for all m = 0:l and given mm = -l:l for 
 * argument beta using the recursion method given in Kostelec and
 * Rockmore (2010) (see equations (4.5)-(4.9)).  For l>abs(mm), require the
 * dl plane to be computed already with values for l-1 *and* l-2.
 * Also takes a table of precomputed square roots of integers and
 * signs to avoid recomputing them.
 *
 * \param[in,out] dlm1p1_line Wigner line.  On input this should be initialised
 * to the line computed for el-2.  On output this will be replaced
 * with the computed line for el.
 * \param[in] dl Wigner plane already computed for el-1.  
 * \param[in] beta Angle to compute Wigner line for.
 * \param[in] L Harmonic band-limit.
 * \param[in] mm Azimuthal harmonic index to compute Wigner line for.
 * \param[in] el Harmonic index to compute Wigner plane for.
 * \param[in] sqrt_tbl Precomputed array of square roots.  The table
 * element at index i should contain the value sqrt(i).  Values from 0
 * to 2*el must be precomputed (i.e. sqrt_tbl should contian 2*el+1
 * elements).
 * \param[in] signs Precomputed array of signs. The array element at
 * index i should contain the value (-1)^i.  Values from 0
 * to el must be precomputed (i.e. signs should contian el+1
 * elements).
 * \retval none
 *
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_beta_kostelec_halfline_table(double *dlm1p1_line, double *dl_line, 
					  double beta, int L, int mm, int el, 
					  double *sqrt_tbl, double *signs) {

  int offset;
  double cosb, sinb, coshb, sinhb;
  double lnAlm, lnAlmm, lnfact2el;
  int m, elm1;
  double elr, elm1r;
  
  // Compute m offset for accessing dl line.
  offset = 0;

  // Compute Wigner plane.
  if (el < abs(mm)) {
    // Do nothing (dl line should remain zero).
    return;
  }
  else if (el == 1) {

    cosb = cos(beta);
    sinb = sin(beta);
    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    if (mm == -1) {
      dlm1p1_line[ 0 + offset] = -sinb / SSHT_SQRT2;
      dlm1p1_line[ 1 + offset] = sinhb * sinhb;
    }
    else if (mm == 0) {
      dlm1p1_line[ 0 + offset] = cosb;
      dlm1p1_line[ 1 + offset] = -sinb / SSHT_SQRT2;
    }
    else {
      // mm == +1
      dlm1p1_line[ 0 + offset] = sinb / SSHT_SQRT2;
      dlm1p1_line[ 1 + offset] = coshb * coshb;
    }

  }
  else if (el == abs(mm)) {

    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    // Initalise line using equation (4.8) or (4.9) from K&R (2010).
    if (mm >= 0) {

      // Initialise using equation (4.8), i.e. top line.
      lnfact2el = logfact(2*el);
      for (m=0; m<=el; m++) {
	lnAlm = (lnfact2el - logfact(el+m) - logfact(el-m)) / 2.0;
	dlm1p1_line[m + offset] = 
	  exp(lnAlm + (el+m)*log(coshb) + (el-m)*log(sinhb));
      }

    }
    else {

      // Initialise using equation (4.9), i.e. bottom line.
      lnfact2el = logfact(2*el);
      for (m=0; m<=el; m++) {
	lnAlm = (lnfact2el - logfact(el+m) - logfact(el-m)) / 2.0;
	dlm1p1_line[m + offset] = 
	  signs[el] * signs[abs(m)]
	  * exp(lnAlm + (el-m)*log(coshb) + (el+m)*log(sinhb));
      }

    }

  }
  else {

    cosb = cos(beta);
    coshb = cos(beta / 2.0);
    sinhb = sin(beta / 2.0);

    elr = (double) el;
    elm1 = el - 1;
    elm1r = (double) elm1;

    // Recuse over line.
    for (m=0; m<=el-1; m++) {

      // Compute 3-term recursion.
      dlm1p1_line[m + offset] = 
	(cosb - m*mm/(elm1r*elr)) * dl_line[m + offset] 
	- 
	sqrt_tbl[elm1+m] * sqrt_tbl[elm1-m] * sqrt_tbl[elm1+mm] * sqrt_tbl[elm1-mm] 
	/ (elm1r * (2.0*elm1r + 1.0))
	* dlm1p1_line[m + offset];
      
      // Perform scaling.
      dlm1p1_line[m + offset] *= 
	el * (2*elm1 + 1.0)
	/ (sqrt_tbl[el-m] * sqrt_tbl[el+m] * sqrt_tbl[el-mm] * sqrt_tbl[el+mm]);

    }

    // Compute edges...
    lnfact2el = logfact(2*el);
    lnAlmm = (lnfact2el - logfact(el+mm) - logfact(el-mm)) / 2.0;
       
    // Right edge.
    dlm1p1_line[el + offset] =  
      signs[el] * signs[abs(mm)]
      * exp(lnAlmm + (el+mm)*log(coshb) + (el-mm)*log(sinhb));

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
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_halfpi_trapani_eighth_table(double *dl, int L, 
					 ssht_dl_size_t dl_size,
					 int el, double *sqrt_tbl) {

  int m, mm, offset, stride;
  double *dmm;
  double t1, t2, s1, s2;

  // Allocate temporary memory.
  dmm = (double*)calloc(el+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dmm)

  // Get mm offset and stride for accessing dl data.
  offset = ssht_dl_get_offset(L, dl_size);
  stride = ssht_dl_get_stride(L, dl_size);

  // Compute Wigner plane.
  if (el == 0) {
    
    dl[0*stride + 0 + offset] = 1.0;

  }
  else {

    // Eqn (9) of T&N (2006).
    dmm[0] = - sqrt_tbl[2*el-1] / sqrt_tbl[2*el]
      * dl[(el-1)*stride + 0 + offset];

    // Eqn (10) of T&N (2006).
    for (mm=1; mm<=el; mm++) {
      dmm[mm] = sqrt_tbl[el] / SSHT_SQRT2 
        * sqrt_tbl[2*el-1] / sqrt_tbl[el+mm] / sqrt_tbl[el+mm-1]
	* dl[(el-1)*stride + (mm-1) + offset];
    }

    // Initialise dl for next el.
    for (mm=0; mm<=el; mm++) {     
      dl[el*stride + mm + offset] = dmm[mm];
    }

/*  LOGICAL BUT *NOT* MOST EFFICIENT ALGORITHM
    // Eqn (11) of T&N (2006).
    for (mm=0; mm<=el; mm++) {

      // m = el-1 case (t2 = 0).
      m = el-1;
      dl[m*stride + mm + offset] = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	* dl[(m+1)*stride + mm + offset];

      // Remaining m cases.
      for (m=el-2; m>=mm; m--) {
    	t1 = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	  * dl[(m+1)*stride + mm + offset];
    	t2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	  * dl[(m+2)*stride + mm + offset];
    	dl[m*stride + mm + offset] = t1 - t2;
      }

    }
*/

    // Eqn (11) of T&N (2006).
    // OPTIMISED FOR MEMORY ACCESS.
    m = el-1;
    s1 = sqrt_tbl[el-m] * sqrt_tbl[el+m+1];
    for (mm=0; mm<=el; mm++) {
      // m = el-1 case (t2 = 0).
      dl[m*stride + mm + offset] = 2e0 * mm / s1
    	* dl[(m+1)*stride + mm + offset];
    }

    for (m=el-2; m>=0; m--) {
      s1 = sqrt_tbl[el-m] * sqrt_tbl[el+m+1];
      s2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1];
      for (mm=0; mm<=m; mm++) {           
    	t1 = 2e0 * mm / s1
    	  * dl[(m+1)*stride + mm + offset];
    	t2 = s2
    	  * dl[(m+2)*stride + mm + offset];
    	dl[m*stride + mm + offset] = t1 - t2;
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
 * \warning THIS GOES UNSTABLE SOMEWHERE ABOVE A BAND-LIMIT OF 1024.
 * Best to compute eighth of plane and then extend to quarter using
 * symmetry (see ssht_dl_halfpi_trapani_eighth_table and
 * ssht_dl_halfpi_trapani_fill_eighth2quarter_table).
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
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_halfpi_trapani_quarter_table(double *dl, int L, 
					 ssht_dl_size_t dl_size,
					 int el, double *sqrt_tbl) {

  int m, mm, offset, stride;
  double *dmm;
  double t1, t2, s1, s2;

  // Allocate temporary memory.
  dmm = (double*)calloc(el+1, sizeof(double));
  SSHT_ERROR_MEM_ALLOC_CHECK(dmm)

  // Get mm offset and stride for accessing dl data.
  offset = ssht_dl_get_offset(L, dl_size);
  stride = ssht_dl_get_stride(L, dl_size);

  // Compute Wigner plane.
  if (el == 0) {
    
    dl[0*stride + 0 + offset] = 1.0;

  }
  else {

    // Eqn (9) of T&N (2006).
    dmm[0] = - sqrt_tbl[2*el-1] / sqrt_tbl[2*el]
      * dl[(el-1)*stride + 0 + offset];

    // Eqn (10) of T&N (2006).
    for (mm=1; mm<=el; mm++) {
      dmm[mm] = sqrt_tbl[el] / SSHT_SQRT2 
        * sqrt_tbl[2*el-1] / sqrt_tbl[el+mm] / sqrt_tbl[el+mm-1]
	* dl[(el-1)*stride + (mm-1) + offset];
    }

    // Initialise dl for next el.
    for (mm=0; mm<=el; mm++) {     
      dl[el*stride + mm + offset] = dmm[mm];
    }

/*  LOGICAL BUT *NOT* MOST EFFICIENT ALGORITHM
    // Eqn (11) of T&N (2006).
    for (mm=0; mm<=el; mm++) {

      // m = el-1 case (t2 = 0).
      m = el-1;
      dl[m*stride + mm + offset] = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	* dl[(m+1)*stride + mm + offset];

      // Remaining m cases.
      for (m=el-2; m>=mm; m--) {
    	t1 = 2e0 * mm / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	  * dl[(m+1)*stride + mm + offset];
    	t2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1]
    	  * dl[(m+2)*stride + mm + offset];
    	dl[m*stride + mm + offset] = t1 - t2;
      }

    }
*/

    // Eqn (11) of T&N (2006).
    // OPTIMISED FOR MEMORY ACCESS.
    m = el-1;
    s1 = sqrt_tbl[el-m] * sqrt_tbl[el+m+1];
    for (mm=0; mm<=el; mm++) {
      // m = el-1 case (t2 = 0).
      dl[m*stride + mm + offset] = 2e0 * mm / s1
    	* dl[(m+1)*stride + mm + offset];
    }

    for (m=el-2; m>=0; m--) {
      s1 = sqrt_tbl[el-m] * sqrt_tbl[el+m+1];
      s2 = sqrt_tbl[el-m-1] * sqrt_tbl[el+m+2] / sqrt_tbl[el-m] / sqrt_tbl[el+m+1];
      for (mm=0; mm<=el; mm++) {
    	t1 = 2e0 * mm / s1
    	  * dl[(m+1)*stride + mm + offset];
    	t2 = s2
    	  * dl[(m+2)*stride + mm + offset];
    	dl[m*stride + mm + offset] = t1 - t2;
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
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_halfpi_trapani_fill_eighth2righthalf_table(double *dl, int L,
							ssht_dl_size_t dl_size,
							int el, double *signs) {

  int m, mm, offset, stride;

  // Get mm offset and stride for accessing dl data.
  offset = ssht_dl_get_offset(L, dl_size);
  stride = ssht_dl_get_stride(L, dl_size);

  // Diagonal symmetry to fill in quarter.
  for (m=0; m<=el; m++)
    for (mm=m+1; mm<=el; mm++)
      dl[m*stride + mm + offset] =
  	signs[m] * signs[mm] * dl[mm*stride + m + offset];

  // Symmetry in mm to fill in half.
  for (m=0; m<=el; m++)
    for (mm=-el; mm<=-1; mm++)
      dl[m*stride + mm + offset] = 
	signs[el] * signs[m] * dl[m*stride - mm + offset];

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
 * \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
 */
void ssht_dl_halfpi_trapani_fill_eighth2quarter_table(double *dl, int L,
						     ssht_dl_size_t dl_size,
						     int el, double *signs) {

  int m, mm, offset, stride;

  // Get mm offset and stride for accessing dl data.
  offset = ssht_dl_get_offset(L, dl_size);
  stride = ssht_dl_get_stride(L, dl_size);

  // Diagonal symmetry to fill in quarter.
  for (m=0; m<=el; m++)
    for (mm=m+1; mm<=el; mm++)
      dl[m*stride + mm + offset] = 
	signs[m] * signs[mm] * dl[mm*stride + m + offset];

}


/*!  
 * Computes the natural logarithm of an (integer) factorial.
 *
 * \param[in] n Integer to compute factorial of.
 * \retval logfactn Natural logarithm of factorial value computed.
 *
 * \author Numerical recipes.
 */
double logfact(int n) {
  
  double y, temp, sum, c[6], loggamma, x;
  int nn;

  if (n < 0) {

    SSHT_ERROR_GENERIC("Factorial argument negative") 

  }
  else {

    // The engine of this function actually calculates the gamma function,
    // for which the real argument is x = n + 1.

    x = (double) (n)  + 1.0;

    // Table of fitting constants.

    c[0] = 76.18009172947146;
    c[1] = - 86.50532032941677;
    c[2] = 24.01409824083091;
    c[3] = - 1.231739572450155;
    c[4] = 0.1208650973866179e-2;
    c[5] = - 0.5395239384953e-5;

    // Add up fit.

    temp = x + 5.5 - (x + 0.5) * log(x + 5.5);
    sum = 1.000000000190015;
    y = x;

    for (nn=0; nn<=5; nn++) {
      y = y + 1.0;
      sum = sum + c[nn] / y;
    }

    loggamma = - temp + log(2.5066282746310005 * sum / x);

  }

  // Finally make explicit the conversion back to log of the factorial.
  return loggamma;

}
