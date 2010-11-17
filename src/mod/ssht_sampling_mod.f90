!------------------------------------------------------------------------------
! ssht_sampling_mod  -- SSHT library sampling class
! 
!! Functionality to define sample positions for various algorithms and to 
!! convert 1D and 2D harmonic indices.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 November 2010
!
! Revisions:
!   November 2010 - Written by Jason McEwen
!------------------------------------------------------------------------------

module ssht_sampling_mod

  use ssht_types_mod
  use ssht_error_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
       ssht_sampling_elm2ind, &
       ssht_sampling_ind2elm, &
       ssht_sampling_dh_t2theta, &
       ssht_sampling_dh_p2phi, &
       ssht_sampling_mw_t2theta, &
       ssht_sampling_mw_p2phi, &
       ssht_sampling_mweo_t2theta, &
       ssht_sampling_mweo_p2phi


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! Flag to indicate Driscoll and Healy sampling.
  integer, public, parameter :: SSHT_SAMPLING_DH = 1

  !! Flag to indicate McEwen and Wiaux sampling.
  integer, public, parameter :: SSHT_SAMPLING_MW = 2

  !! Flag to indicate McEwen and Wiaux sampling for even-odd algorithm.
  integer, public, parameter :: SSHT_SAMPLING_MWEO = 3

  !! Flag to indicate default sampling (McEwen and Wiaux sampling).
  integer, public, parameter :: SSHT_SAMPLING_DEFAULT = SSHT_SAMPLING_MW


  !---------------------------------------
  ! Data types
  !---------------------------------------

  ! None.


  !----------------------------------------------------------------------------

contains


  !============================================================================
  ! Sampling relations
  !============================================================================


  !----------------------------------------------------------------------------
  ! ssht_sampling_dh_t2theta
  !
  !! Covert theta index to angle for Driscoll and Healy sampling.
  !!
  !! Notes:
  !!  - t ranges from [0 .. 2*L-1] => 2*L points in (0,pi).
  !!
  !! Variables:
  !!  - t: Theta index [input].
  !!  - theta: Theta angle [output].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------
  
  function ssht_sampling_dh_t2theta(t, L) result (theta)

    integer, intent(in) :: t
    integer, intent(in) :: L
    real(dp) :: theta

    theta = (2d0*t+1d0)*PI / (4d0*L)

  end function ssht_sampling_dh_t2theta


  !----------------------------------------------------------------------------
  ! ssht_sampling_dh_p2phi
  !
  !! Covert phi index to angle for Driscoll and Healy sampling.
  !!
  !! Notes:
  !!  - p ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
  !!
  !! Variables:
  !!  - p: Phi index [input].
  !!  - phi: Phi angle [output].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------
  
  function ssht_sampling_dh_p2phi(p, L) result (phi)

    integer, intent(in) :: p
    integer, intent(in) :: L
    real(dp) :: phi

    phi = 2d0*p*PI / (2d0*L - 1d0)	

  end function ssht_sampling_dh_p2phi


  !----------------------------------------------------------------------------
  ! ssht_sampling_mw_t2theta
  !
  !! Covert theta index to angle for McEwen and Wiaux sampling.
  !!
  !! Notes:
  !!  - t ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
  !!
  !! Variables:
  !!  - t: Theta index [input].
  !!  - theta: Theta angle [output].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------
 
  function ssht_sampling_mw_t2theta(t, L) result (theta)

    integer, intent(in) :: t
    integer, intent(in) :: L
    real(dp) :: theta

    theta = (2d0*t+1d0)*PI / (2d0*L - 1d0)

  end function ssht_sampling_mw_t2theta


  !----------------------------------------------------------------------------
  ! ssht_sampling_mw_p2phi
  !
  !! Covert phi index to angle for McEwen and Wiaux sampling.
  !!
  !! Notes:
  !!  - p ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
  !!
  !! Variables:
  !!  - p: Phi index [input].
  !!  - phi: Phi angle [output].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  function ssht_sampling_mw_p2phi(p, L) result (phi)

    integer, intent(in) :: p
    integer, intent(in) :: L
    real(dp) :: phi

    phi = 2d0*p*PI / (2d0*L - 1d0)	

  end function ssht_sampling_mw_p2phi


  !----------------------------------------------------------------------------
  ! ssht_sampling_mweo_t2theta
  !
  !! Covert theta index to angle for McEwen and Wiaux even-odd sampling.
  !!
  !! Notes:
  !!  - t ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
  !!
  !! Variables:
  !!  - t: Theta index [input].
  !!  - theta: Theta angle [output].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  function ssht_sampling_mweo_t2theta(t, L) result (theta)

    integer, intent(in) :: t
    integer, intent(in) :: L
    real(dp) :: theta

    theta = (2d0*t+1d0)*PI / (2d0*L - 1d0)

  end function ssht_sampling_mweo_t2theta


  !----------------------------------------------------------------------------
  ! ssht_sampling_mweo_p2phi
  !
  !! Covert phi index to angle for McEwen and Wiaux even-odd sampling.
  !!
  !! Notes:
  !!  - p ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
  !!
  !! Variables:
  !!  - p: Phi index [input].
  !!  - phi: Phi angle [output].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  function ssht_sampling_mweo_p2phi(p, L) result (phi)

    integer, intent(in) :: p
    integer, intent(in) :: L
    real(dp) :: phi

    phi = (2d0*p+1d0)*PI / (2d0*L - 1d0)

  end function ssht_sampling_mweo_p2phi


  !============================================================================
  ! Harmonic index relations
  !============================================================================

  
  !----------------------------------------------------------------------------
  ! ssht_sampling_elm2ind
  !
  !! Covert (el,m) harmonic indices to 1D index used to access flm array.
  !!
  !! Notes:
  !!  - el ranges from [0 .. L-1].
  !!  - m ranges from [-el .. el].
  !!  - ind ranges from [0 .. L**2-1]
  !!
  !! Variables:
  !!  - ind: 1D index to access flm array [output].
  !!  - el: Harmonic index [input].
  !!  - m: Azimuthal harmonic index [input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_sampling_elm2ind(ind, el, m)

    integer, intent(out) :: ind
    integer, intent(in) :: el, m

    ind = el**2 + el + m

  end subroutine ssht_sampling_elm2ind


  !----------------------------------------------------------------------------
  ! ssht_sampling_ind2elm
  !
  !! Covert 1D index used to access flm array to (el,m) harmonic indices.
  !!
  !! Notes:
  !!  - el ranges from [0 .. L-1].
  !!  - m ranges from [-el .. el].
  !!  - ind ranges from [0 .. L**2-1]
  !!
  !! Variables:
  !!  - el: Harmonic index [output].
  !!  - m: Azimuthal harmonic index [output].
  !!  - ind: 1D index to access flm array [input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_sampling_ind2elm(el, m, ind)

    integer, intent(out) :: el, m
    integer, intent(in) :: ind

    el = floor(sqrt(real(ind,dp)))
    m = ind - el**2 - el

  end subroutine ssht_sampling_ind2elm


end module ssht_sampling_mod
