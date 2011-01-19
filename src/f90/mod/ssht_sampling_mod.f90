!------------------------------------------------------------------------------
! ssht_sampling_mod  -- SSHT library sampling class
! 
!! Functionality to define sample positions for various algorithms, to compute
!! weights and to convert 1D and 2D harmonic indices.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
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
       ssht_sampling_gl_p2phi, &
       ssht_sampling_mw_t2theta, &
       ssht_sampling_mw_p2phi, &
       ssht_sampling_mweo_t2theta, &
       ssht_sampling_mweo_p2phi, &
       ssht_sampling_weight_dh, &
       ssht_sampling_weight_mw, &
       ssht_sampling_gl_thetas_weights


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  !! Flag to indicate Driscoll and Healy sampling.
  integer, public, parameter :: SSHT_SAMPLING_DH = 1

  !! Flag to indicate Gauss-Legendre sampling.
  integer, public, parameter :: SSHT_SAMPLING_GL = 2

  !! Flag to indicate McEwen and Wiaux sampling.
  integer, public, parameter :: SSHT_SAMPLING_MW = 3

  !! Flag to indicate McEwen and Wiaux sampling for even-odd algorithm.
  integer, public, parameter :: SSHT_SAMPLING_MWEO = 4

  !! Flag to indicate default sampling (McEwen and Wiaux sampling).
  integer, public, parameter :: SSHT_SAMPLING_DEFAULT = SSHT_SAMPLING_MW


  !---------------------------------------
  ! Data types
  !---------------------------------------

  ! None.


  !----------------------------------------------------------------------------

contains


  !============================================================================
  ! Sampling weights
  !============================================================================


  !--------------------------------------------------------------------------
  ! ssht_sampling_weight_dh
  !
  !! Compute Discoll and Healy weights.
  !!
  !! Variables:
  !!  - theta_t: Theta value to compute weight for [input].
  !!  - L: Harmonic band-limit [input].
  !!  - w: Corresponding weight [output]
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2007
  !
  ! Revisions:
  !   October 2007 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function ssht_sampling_weight_dh(theta_t, L) result(w)

    real(dp), intent(in) :: theta_t
    integer, intent(in) :: L
    real(dp) :: w	

    integer :: k

    w = 0d0
    do k = 0,L-1
       w = w + sin((2d0*k+1d0)*theta_t) / real(2d0*k+1d0,dp)
    end do
    w = (2d0/real(L,dp)) * sin(theta_t) * w

  end function ssht_sampling_weight_dh


  !--------------------------------------------------------------------------
  ! ssht_sampling_weight_mw
  !
  !! Compute weights for toroidal extension.
  !!
  !! Variables:
  !!  - p: Integer index to compute weight for [input].
  !!  - w: Corresponding weight [output]
  !
  !! @author J. D. McEwen
  !! @version 0.1 October 2010
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  function ssht_sampling_weight_mw(p) result(w)

    integer, intent(in) :: p
    complex(dpc) :: w

    if(p == 1) then
       w = I * PION2
    elseif(p == -1) then
       w = - I * PION2
    elseif(mod(p,2) == 0 ) then
       ! Even case
       w = 2d0 / (1d0 - p**2)
    else
       ! Odd case (|p| /= 1)
       w = 0d0
    end if

  end function ssht_sampling_weight_mw


  !--------------------------------------------------------------------------
  ! ssht_sampling_gl_thetas_weights
  !
  !! Compute theta positions (roots of Legendre polynomials) and 
  !! corresponding weights.
  !!
  !! Variables:
  !!  - thetas(0:L-1): Theta positions [output].
  !!  - weights(0:L-1): Corresponding weights [output]
  !!  - L: Harmonic band-limit [input].
  !
  !! @author J. D. McEwen
  !! @version 0.1 November 2010
  !
  ! Revisions:
  !   November 2010 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine ssht_sampling_gl_thetas_weights(thetas, weights, L)

    integer, intent(in) :: L
    real(dp), intent(out) :: thetas(0:L-1)
    real(dp), intent(out) :: weights(0:L-1)
    
    real(dp) :: znodes(0:L-1)
    integer :: p

    call gauleg(-1d0, 1d0, znodes(0:L-1), weights(0:L-1), L)

    do p = 0, L-1
       thetas(p) = acos(znodes(p))
    end do

  end subroutine ssht_sampling_gl_thetas_weights


  !--------------------------------------------------------------------------
  ! gauleg
  !
  !! Given the lower and upper limits of integration x1 and x2, this
  !! routine returns arrays x[1..n] and w[1..n] of length n,
  !! containing the abscissas and weights of the Gauss-Legendre
  !! n-point quadrature formula.
  !!
  !! Variables:
  !!   - x1: Lower bound of range [input].
  !!   - x2: Upper bound of range [input].
  !!   - x: Node positions (i.e. roots of Legendre polynomials) [output].
  !!   - w: Corresponding weights [output].
  !!   - n: Number of points [input].
  !!   - integral: Value of the evaluated integral [output].
  !
  !! @author Numerical recipes.
  !
  ! Revisions:
  !   November 2010 - Adapted from Numerical Recipes by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine gauleg(x1, x2, x, w, n)

    integer, intent(in) :: n
    real(dp), intent(in) :: x1,x2
    real(dp), intent(out) :: x(1:n),w(1:n)

    real(dp), parameter :: EPS = 1d-14
    integer :: i,j,m
    real(dp) :: p1,p2,p3,pp,xl,xm,z,z1

    m=(n+1)/2
    xm=0.5d0*(x2+x1)
    xl=0.5d0*(x2-x1)
    do i=1,m
       z=cos(3.141592654d0*(i-.25d0)/(n+.5d0))
       do 
          p1=1.d0
          p2=0.d0
          do j=1,n
             p3=p2
             p2=p1
             p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/j
          end do
          pp=n*(z*p1-p2)/(z*z-1.d0)
          z1=z
          z=z1-p1/pp
          if(abs(z-z1) <= EPS) exit
       end do
       x(i)=xm-xl*z
       x(n+1-i)=xm+xl*z
       w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
       w(n+1-i)=w(i)
    end do

  end subroutine gauleg


  !============================================================================
  ! Sampling relations
  !============================================================================


  !----------------------------------------------------------------------------
  ! ssht_sampling_dh_t2theta
  !
  !! Convert theta index to angle for Driscoll and Healy sampling.
  !!
  !! Notes:
  !!  - t ranges from [0 .. 2*L-1] => 2*L points in (0,pi).
  !!
  !! Variables:
  !!  - t: Theta index [input].
  !!  - L: Harmonic band-limit [input].
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
  !! Convert phi index to angle for Driscoll and Healy sampling.
  !!
  !! Notes:
  !!  - p ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
  !!
  !! Variables:
  !!  - p: Phi index [input].
  !!  - L: Harmonic band-limit [input].
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
  ! ssht_sampling_gl_p2phi
  !
  !! Convert phi index to angle for Gauss-Legendre sampling.
  !!
  !! Notes:
  !!  - p ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
  !!
  !! Variables:
  !!  - p: Phi index [input].
  !!  - L: Harmonic band-limit [input].
  !!  - phi: Phi angle [output].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   November 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------
  
  function ssht_sampling_gl_p2phi(p, L) result (phi)

    integer, intent(in) :: p
    integer, intent(in) :: L
    real(dp) :: phi

    phi = 2d0*p*PI / (2d0*L - 1d0)	

  end function ssht_sampling_gl_p2phi


  !----------------------------------------------------------------------------
  ! ssht_sampling_mw_t2theta
  !
  !! Convert theta index to angle for McEwen and Wiaux sampling.
  !!
  !! Notes:
  !!  - t ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
  !!
  !! Variables:
  !!  - t: Theta index [input].
  !!  - L: Harmonic band-limit [input].
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
  !! Convert phi index to angle for McEwen and Wiaux sampling.
  !!
  !! Notes:
  !!  - p ranges from [0 .. 2*L-2] => 2*L-1 points in [0,2*pi).
  !!
  !! Variables:
  !!  - p: Phi index [input].
  !!  - L: Harmonic band-limit [input].
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
  !! Convert theta index to angle for McEwen and Wiaux even-odd sampling.
  !!
  !! Notes:
  !!  - t ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
  !!
  !! Variables:
  !!  - t: Theta index [input].
  !!  - L: Harmonic band-limit [input].
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
  !! Convert phi index to angle for McEwen and Wiaux even-odd sampling.
  !!
  !! Notes:
  !!  - p ranges from [0 .. 2*L-2] => 2*L-1 points in (0,2*pi).
  !!
  !! Variables:
  !!  - p: Phi index [input].
  !!  - L: Harmonic band-limit [input].
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
  !! Convert (el,m) harmonic indices to 1D index used to access flm array.
  !!
  !! Notes:
  !!  - el ranges from [0 .. L-1].
  !!  - m ranges from [-el .. el].
  !!  - ind ranges from [0 .. L**2-1].
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
  !! Convert 1D index used to access flm array to (el,m) harmonic indices.
  !!
  !! Notes:
  !!  - el ranges from [0 .. L-1].
  !!  - m ranges from [-el .. el].
  !!  - ind ranges from [0 .. L**2-1].
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
