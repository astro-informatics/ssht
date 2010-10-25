!------------------------------------------------------------------------------
! ssht_types_mod -- SSHT library types class
!
!! Definition of intrinsic types and constants used in the ssht library.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 October 2007
!
! Revisions:
!   October 2007 - Written by Jason McEwen
!------------------------------------------------------------------------------

module ssht_types_mod

  implicit none

  private


  ! --------------------------------------
  ! Intrinsic type definitions
  ! --------------------------------------

  !! Type definition for single precision real.
  integer, public, parameter :: sp  = SELECTED_REAL_KIND(5,30)

  !! Type definition for double precision real.
  integer, public, parameter :: dp  = SELECTED_REAL_KIND(12,200)

  !! Type definition for single precisison complex.
  integer, public, parameter :: spc = KIND((1.0_sp, 1.0_sp))

  !! Type definition for double precision complex.
  integer, public, parameter :: dpc = KIND((1.0_dp, 1.0_dp))


  ! --------------------------------------
  ! Tolerances
  ! --------------------------------------

	!! Admissibility tolerance.
	real(dp), public, parameter :: TOL_ADMISS = 1e-12

	!! Numerical integration tolerance.
	real(dp), public, parameter :: TOL_QUAD = 1d-12

	!! Limit opening tolerance for definite numerical integration.
	real(dp), public, parameter :: TOL_LIMIT = 1d-12

	!! Tolerance for floor and ceiling functions.
	real(dp), public, parameter :: TOL_CEIL = 1d-5
	


  ! --------------------------------------
  ! Constants
  ! --------------------------------------

  !! String buffer length.
  integer, public, parameter :: STRING_LEN = 256

  !! PI definition.
  real(dp), public, parameter :: PI = 3.141592653589793238462643383279502884197_dp

  !! PI/2 definition.
  real(dp), public, parameter :: PION2 = 1.570796326794896619231321691639751442099_dp

  !! Complex unit definition.
  complex(dpc), public, parameter :: I = (0.0_dp, 1.0_dp)


end module ssht_types_mod
