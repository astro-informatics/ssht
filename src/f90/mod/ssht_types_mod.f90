! SSHT package to perform spin spherical harmonic transforms
! Copyright (C) 2011  Jason McEwen
! See LICENSE.txt for license details


!> \mainpage SSHT Fortran documentation
!!
!! The SSHT code provides functionality to perform fast and exact spin
!! spherical harmonic transforms based on the sampling theorem on the
!! sphere derived in our paper: <i>A novel sampling theorem on the
!! sphere</i> (<a href="http://arxiv.org/abs/XXX.XXX">ArXiv</a>|
!! <a href="http://dx.doi.org/10.1111/XXX">DOI</a>).
!!   
!! <p>
!! We document the Fortran source code here.  For an example of usage, 
!! see the ssht_test.f90 program.
!! For installation instructions, see the general SSHT 
!! documentation available 
!! <a href="../../index_ssht.html">here</a>.
!!
!! \author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!!

!------------------------------------------------------------------------------
! ssht_types_mod -- SSHT library types class
!
!> Definition of intrinsic types and constants used in the ssht library.
!!
!! @author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!
! Revisions:
!   October 2010 - Written by Jason McEwen
!------------------------------------------------------------------------------

module ssht_types_mod

  implicit none

  private


  ! --------------------------------------
  ! Intrinsic type definitions
  ! --------------------------------------

  !> Type definition for single precision real.
  integer, public, parameter :: sp  = SELECTED_REAL_KIND(5,30)

  !> Type definition for double precision real.
  integer, public, parameter :: dp  = SELECTED_REAL_KIND(12,200)

  !> Type definition for single precisison complex.
  integer, public, parameter :: spc = KIND((1.0_sp, 1.0_sp))

  !> Type definition for double precision complex.
  integer, public, parameter :: dpc = KIND((1.0_dp, 1.0_dp))


  ! --------------------------------------
  ! Tolerances
  ! --------------------------------------

  !> Realitive reality tolerance.
  real(dp), public, parameter :: TOL_RELATIVE_REALITY = 1e-6


  ! --------------------------------------
  ! Constants
  ! --------------------------------------

  !> String buffer length.
  integer, public, parameter :: STRING_LEN = 256

  !> PI definition.
  real(dp), public, parameter :: PI = 3.141592653589793238462643383279502884197_dp

  !> PI/2 definition.
  real(dp), public, parameter :: PION2 = 1.570796326794896619231321691639751442099_dp

  !> SQRT(2) definition.
  real(dp), public, parameter :: SQRT2 = 1.41421356237309504880168872420969807856967_dp

  !> Complex unit definition.
  complex(dpc), public, parameter :: I = (0.0_dp, 1.0_dp)


  ! --------------------------------------
  ! Prompt
  ! --------------------------------------

  !> Verbosity prompt.
  character(len=*), public, parameter :: SSHT_PROMPT = '[ssht] '


end module ssht_types_mod
