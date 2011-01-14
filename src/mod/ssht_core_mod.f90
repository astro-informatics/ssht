!------------------------------------------------------------------------------
! ssht_core_mod  -- SSHT library core class
! 
!! Functionality to perform spin spherical harmonic transform on the sphere.
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   October 2010 - Written by Jason McEwen
!------------------------------------------------------------------------------

module ssht_core_mod

  use ssht_types_mod
  use ssht_error_mod
  use ssht_dl_mod
  use ssht_sampling_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
       ssht_core_dh_inverse, &
       ssht_core_dh_forward, &
       ssht_core_gl_inverse, &
       ssht_core_gl_forward, &
       ssht_core_mweo_inverse, &
       ssht_core_mweo_forward, &
       ssht_core_mw_inverse, &
       ssht_core_mw_forward, &
       ssht_core_dh_inverse_real, &
       ssht_core_dh_forward_real, &
       ssht_core_gl_inverse_real, &
       ssht_core_gl_forward_real, &
       ssht_core_mweo_inverse_real, &
       ssht_core_mweo_forward_real, &
       ssht_core_mw_inverse_real, &
       ssht_core_mw_forward_real


  !----------------------------------------------------------------------------
  ! DH interfaces
  !----------------------------------------------------------------------------

  ! Define implementations available publicly.

  interface ssht_core_dh_inverse
     module procedure ssht_core_dh_inverse_sov_sym
  end interface

  interface ssht_core_dh_forward
     module procedure ssht_core_dh_forward_sov_sym
  end interface

  interface ssht_core_dh_inverse_real
     module procedure ssht_core_dh_inverse_sov_sym_real
  end interface

  interface ssht_core_dh_forward_real
     module procedure ssht_core_dh_forward_sov_sym_real
  end interface

  ! Define implementations to use in south pole wrapper routines.

  interface ssht_core_dh_inverse_rsp
     module procedure ssht_core_dh_inverse_sov_sym
  end interface

  interface ssht_core_dh_forward_rsp
     module procedure ssht_core_dh_forward_sov_sym
  end interface

  interface ssht_core_dh_inverse_real_rsp
     module procedure ssht_core_dh_inverse_sov_sym_real
  end interface

  interface ssht_core_dh_forward_real_rsp
     module procedure ssht_core_dh_forward_sov_sym_real
  end interface


  !----------------------------------------------------------------------------
  ! GL interfaces
  !----------------------------------------------------------------------------

  ! Define implementations available publicly.

  interface ssht_core_gl_inverse
     module procedure ssht_core_gl_inverse_sov_sym
  end interface

  interface ssht_core_gl_forward
     module procedure ssht_core_gl_forward_sov_sym
  end interface

  interface ssht_core_gl_inverse_real
     module procedure ssht_core_gl_inverse_sov_sym_real
  end interface

  interface ssht_core_gl_forward_real
     module procedure ssht_core_gl_forward_sov_sym_real
  end interface

  ! Define implementations to use in south pole wrapper routines.

  interface ssht_core_gl_inverse_rsp
     module procedure ssht_core_gl_inverse_sov_sym
  end interface

  interface ssht_core_gl_forward_rsp
     module procedure ssht_core_gl_forward_sov_sym
  end interface

  interface ssht_core_gl_inverse_real_rsp
     module procedure ssht_core_gl_inverse_sov_sym_real
  end interface

  interface ssht_core_gl_forward_real_rsp
     module procedure ssht_core_gl_forward_sov_sym_real
  end interface


  !----------------------------------------------------------------------------
  ! MWEO interfaces
  !----------------------------------------------------------------------------

  ! Define implementations available publicly.

  interface ssht_core_mweo_inverse
     module procedure ssht_core_mweo_inverse_sp
  end interface

  interface ssht_core_mweo_forward
     module procedure ssht_core_mweo_forward_sp
  end interface

  interface ssht_core_mweo_inverse_real
     module procedure ssht_core_mweo_inverse_real_sp
  end interface

  interface ssht_core_mweo_forward_real
     module procedure ssht_core_mweo_forward_real_sp
  end interface

  ! Define implementations to use in south pole wrapper routines.

  interface ssht_core_mweo_inverse_rsp
     module procedure ssht_core_mweo_inverse_sov_sym
  end interface

  interface ssht_core_mweo_forward_rsp
     module procedure ssht_core_mweo_forward_sov_conv_sym
  end interface

  interface ssht_core_mweo_inverse_real_rsp
     module procedure ssht_core_mweo_inverse_sov_sym_real
  end interface

  interface ssht_core_mweo_forward_real_rsp
     module procedure ssht_core_mweo_forward_sov_conv_sym_real
  end interface


  !----------------------------------------------------------------------------
  ! MW interfaces
  !----------------------------------------------------------------------------

  ! Define implementations available publicly.

  interface ssht_core_mw_inverse
     module procedure ssht_core_mw_inverse_sp
  end interface

  interface ssht_core_mw_forward
     module procedure ssht_core_mw_forward_sp
  end interface

  interface ssht_core_mw_inverse_real
     module procedure ssht_core_mw_inverse_real_sp
  end interface

  interface ssht_core_mw_forward_real
     module procedure ssht_core_mw_forward_real_sp
  end interface

  ! Define implementations to use in south pole wrapper routines.

  interface ssht_core_mw_inverse_rsp
     module procedure ssht_core_mw_inverse_sov_sym
  end interface

  interface ssht_core_mw_forward_rsp
     module procedure ssht_core_mw_forward_sov_conv_sym
  end interface

  interface ssht_core_mw_inverse_real_rsp
     module procedure ssht_core_mw_inverse_sov_sym_real
  end interface

  interface ssht_core_mw_forward_real_rsp
     module procedure ssht_core_mw_forward_sov_conv_sym_real
  end interface


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  integer, parameter :: FFTW_ESTIMATE=64, FFTW_MEASURE=0
  integer, parameter :: FFTW_FORWARD=-1, FFTW_BACKWARD=1


  !---------------------------------------
  ! Data types
  !---------------------------------------

  ! None.


  !----------------------------------------------------------------------------

contains


  !============================================================================
  ! South pole interface routines
  !============================================================================


  !----------------------------------------------------------------------------
  ! DH
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_dh_inverse_sp
  !
  !! South pole wrapper for inverse transform using DH method.
  !!
  !! Variables:
  !!   - f(0:2*L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - phi_sp: Phi angle corresponding to quoted sample at South pole 
  !!     [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------


  subroutine ssht_core_dh_inverse_sp(f, f_sp, phi_sp, flm, L, spin, verbosity) 

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: phi_sp
    complex(dpc), intent(out) :: f_sp
    complex(dpc), intent(out) :: f(0:2*L-2, 0:2*L-2)
   
    complex(dpc) :: f_ext(0:2*L-1, 0:2*L-2)

    call ssht_core_dh_inverse_rsp(f_ext, flm, L, spin, verbosity)

    f(0:2*L-2, 0:2*L-2) = f_ext(0:2*L-2, 0:2*L-2)
    f_sp = f_ext(2*L-1, 0)
    phi_sp = ssht_sampling_dh_p2phi(0, L)

  end subroutine ssht_core_dh_inverse_sp


  !----------------------------------------------------------------------------
  ! ssht_core_dh_forward_sp
  !
  !! South pole wrapper for forward transform using DH method.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:2*L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [input].
  !!   - f_sp: Function on sphere South pole sample [input].
  !!   - phi_sp: Phi angle corresponding to quoted sample at South pole 
  !!     [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_forward_sp(flm, f, f_sp, phi_sp, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:2*L-2 ,0:2*L-2)
    real(dp), intent(in) :: phi_sp
    complex(dpc), intent(in) :: f_sp
    complex(dpc), intent(out) :: flm(0:L**2-1)

    complex(dpc) :: f_ext(0:2*L-1, 0:2*L-2)
    integer :: p
    real(dp) :: phi
  
    f_ext(0:2*L-2 ,0:2*L-2) = f(0:2*L-2 ,0:2*L-2)
    do p = 0, 2*L-2
       phi = ssht_sampling_dh_p2phi(p, L)
       f_ext(2*L-1, p) = f_sp * exp(I*spin*(phi-phi_sp))
    end do

    call ssht_core_dh_forward_rsp(flm, f_ext, L, spin, verbosity)

  end subroutine ssht_core_dh_forward_sp


  !----------------------------------------------------------------------------
  ! ssht_core_dh_inverse_real_sp
  !
  !! South pole wrapper for inverse transform of real scalar function
  !! using DH method.
  !!
  !! Variables:
  !!   - f(0:2*L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_inverse_real_sp(f, f_sp, flm, L, verbosity) 

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: f_sp
    real(dp), intent(out) :: f(0:2*L-2, 0:2*L-2)
   
    real(dp) :: f_ext(0:2*L-1, 0:2*L-2)

    call ssht_core_dh_inverse_real_rsp(f_ext, flm, L, verbosity)

    f(0:2*L-2, 0:2*L-2) = f_ext(0:2*L-2, 0:2*L-2)
    f_sp = f_ext(2*L-1, 0)

  end subroutine ssht_core_dh_inverse_real_sp


  !----------------------------------------------------------------------------
  ! ssht_core_dh_forward_real_sp
  !
  !! South pole wrapper for forward transform of real scalar function
  !! using DH method.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - f(0:2*L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - L: Harmonic band-limit [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_forward_real_sp(flm, f, f_sp, L, verbosity)

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    real(dp), intent(in) :: f(0:2*L-2, 0:2*L-2)
    real(dp), intent(in) :: f_sp
    complex(dpc), intent(out) :: flm(0:L**2-1)

    real(dp) :: f_ext(0:2*L-1, 0:2*L-2)

    f_ext(0:2*L-2, 0:2*L-2) = f(0:2*L-2, 0:2*L-2)
    f_ext(2*L-1, 0:2*L-2) = f_sp

    call ssht_core_dh_forward_real_rsp(flm, f_ext, L, verbosity)

  end subroutine ssht_core_dh_forward_real_sp


  !----------------------------------------------------------------------------
  ! GL
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_gl_inverse_sp
  !
  !! South pole wrapper for inverse transform using GL method.
  !!
  !! Variables:
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - phi_sp: Phi angle corresponding to quoted sample at South pole 
  !!     [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_gl_inverse_sp(f, f_sp, phi_sp, flm, L, spin, verbosity) 

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: phi_sp
    complex(dpc), intent(out) :: f_sp
    complex(dpc), intent(out) :: f(0:L-2, 0:2*L-2)
   
    complex(dpc) :: f_ext(0:L-1, 0:L-2)

    call ssht_core_gl_inverse_rsp(f_ext, flm, L, spin, verbosity)

    f(0:L-2, 0:2*L-2) = f_ext(0:L-2, 0:2*L-2)
    f_sp = f_ext(L-1, 0)
    phi_sp = ssht_sampling_gl_p2phi(0, L)

  end subroutine ssht_core_gl_inverse_sp


  !----------------------------------------------------------------------------
  ! ssht_core_gl_forward_sp
  !
  !! South pole wrapper for forward transform using GL method.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [input].
  !!   - f_sp: Function on sphere South pole sample [input].
  !!   - phi_sp: Phi angle corresponding to quoted sample at South pole 
  !!     [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_gl_forward_sp(flm, f, f_sp, phi_sp, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-2 ,0:2*L-2)
    real(dp), intent(in) :: phi_sp
    complex(dpc), intent(in) :: f_sp
    complex(dpc), intent(out) :: flm(0:L**2-1)

    complex(dpc) :: f_ext(0:L-1, 0:2*L-2)
    integer :: p
    real(dp) :: phi
  
    f_ext(0:L-2 ,0:2*L-2) = f(0:L-2 ,0:2*L-2)
    do p = 0, 2*L-2
       phi = ssht_sampling_gl_p2phi(p, L)
       f_ext(L-1, p) = f_sp * exp(I*spin*(phi-phi_sp))
    end do

    call ssht_core_gl_forward_rsp(flm, f_ext, L, spin, verbosity)

  end subroutine ssht_core_gl_forward_sp


  !----------------------------------------------------------------------------
  ! ssht_core_gl_inverse_real_sp
  !
  !! South pole wrapper for inverse transform of real scalar function
  !! using GL method.
  !!
  !! Variables:
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_gl_inverse_real_sp(f, f_sp, flm, L, verbosity) 

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: f_sp
    real(dp), intent(out) :: f(0:L-2, 0:2*L-2)
   
    real(dp) :: f_ext(0:L-1, 0:2*L-2)

    call ssht_core_gl_inverse_real_rsp(f_ext, flm, L, verbosity)

    f(0:L-2, 0:2*L-2) = f_ext(0:L-2, 0:2*L-2)
    f_sp = f_ext(L-1, 0)

  end subroutine ssht_core_gl_inverse_real_sp


  !----------------------------------------------------------------------------
  ! ssht_core_gl_forward_real_sp
  !
  !! South pole wrapper for forward transform of real scalar function
  !! using GL method.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - L: Harmonic band-limit [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_gl_forward_real_sp(flm, f, f_sp, L, verbosity)

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    real(dp), intent(in) :: f(0:L-2, 0:2*L-2)
    real(dp), intent(in) :: f_sp
    complex(dpc), intent(out) :: flm(0:L**2-1)

    real(dp) :: f_ext(0:L-1, 0:2*L-2)

    f_ext(0:L-2, 0:2*L-2) = f(0:L-2, 0:2*L-2)
    f_ext(L-1, 0:2*L-2) = f_sp

    call ssht_core_gl_forward_real_rsp(flm, f_ext, L, verbosity)

  end subroutine ssht_core_gl_forward_real_sp


  !----------------------------------------------------------------------------
  ! MWEO
  !----------------------------------------------------------------------------

 
  !----------------------------------------------------------------------------
  ! ssht_core_mweo_inverse_sp
  !
  !! South pole wrapper for inverse transform using MWEO method.
  !!
  !! Variables:
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - phi_sp: Phi angle corresponding to quoted sample at South pole 
  !!     [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------
  
  subroutine ssht_core_mweo_inverse_sp(f, f_sp, phi_sp, flm, L, spin, verbosity) 

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: phi_sp
    complex(dpc), intent(out) :: f_sp
    complex(dpc), intent(out) :: f(0:L-2, 0:2*L-2)
    
    complex(dpc) :: f_ext(0:L-1, 0:2*L-2)

    call ssht_core_mweo_inverse_rsp(f_ext, flm, L, spin, verbosity)

    f(0:L-2, 0:2*L-2) = f_ext(0:L-2, 0:2*L-2)
    f_sp = f_ext(L-1, 0)
    phi_sp = ssht_sampling_mweo_p2phi(0, L)

  end subroutine ssht_core_mweo_inverse_sp


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_forward_sp
  !
  !! South pole wrapper for forward transform using MWEO method.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [input].
  !!   - f_sp: Function on sphere South pole sample [input].
  !!   - phi_sp: Phi angle corresponding to quoted sample at South pole 
  !!     [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_forward_sp(flm, f, f_sp, phi_sp, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-2 ,0:2*L-2)
    real(dp), intent(in) :: phi_sp
    complex(dpc), intent(in) :: f_sp
    complex(dpc), intent(out) :: flm(0:L**2-1)

    complex(dpc) :: f_ext(0:L-1, 0:2*L-2)
    integer :: p
    real(dp) :: phi
  
    f_ext(0:L-2 ,0:2*L-2) = f(0:L-2 ,0:2*L-2)
    do p = 0, 2*L-2
       phi = ssht_sampling_mweo_p2phi(p, L)
       f_ext(L-1, p) = f_sp * exp(I*spin*(phi-phi_sp))
    end do

    call ssht_core_mweo_forward_rsp(flm, f_ext, L, spin, verbosity)

  end subroutine ssht_core_mweo_forward_sp


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_inverse_real_sp
  !
  !! South pole wrapper for inverse transform of real scalar function
  !! using MWEO method.
  !!
  !! Variables:
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_inverse_real_sp(f, f_sp, flm, L, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: f_sp
    real(dp), intent(out) :: f(0:L-2, 0:2*L-2)

    real(dp) :: f_ext(0:L-1, 0:2*L-2)

    call ssht_core_mweo_inverse_real_rsp(f_ext, flm, L, verbosity)

    f(0:L-2, 0:2*L-2) = f_ext(0:L-2, 0:2*L-2)
    f_sp = f_ext(L-1, 0)

  end subroutine ssht_core_mweo_inverse_real_sp


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_forward_real_sp
  !
  !! South pole wrapper for forward transform of real scalar function
  !! using MWEO method.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - L: Harmonic band-limit [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_forward_real_sp(flm, f, f_sp, L, verbosity)

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    real(dp), intent(in) :: f(0:L-2, 0:2*L-2)
    real(dp), intent(in) :: f_sp
    complex(dpc), intent(out) :: flm(0:L**2-1)

    real(dp) :: f_ext(0:L-1, 0:2*L-2)

    f_ext(0:L-2, 0:2*L-2) = f(0:L-2, 0:2*L-2)
    f_ext(L-1, 0:2*L-2) = f_sp

    call ssht_core_mweo_forward_real_rsp(flm, f_ext, L, verbosity)

  end subroutine ssht_core_mweo_forward_real_sp


  !----------------------------------------------------------------------------
  ! MW
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_mw_inverse_sp
  !
  !! South pole wrapper for inverse transform using MW method.
  !!
  !! Variables:
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - phi_sp: Phi angle corresponding to quoted sample at South pole 
  !!     [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_inverse_sp(f, f_sp, phi_sp, flm, L, spin, verbosity) 

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: phi_sp
    complex(dpc), intent(out) :: f_sp
    complex(dpc), intent(out) :: f(0:L-2, 0:2*L-2)
    
    complex(dpc) :: f_ext(0:L-1, 0:2*L-2)

    call ssht_core_mw_inverse_rsp(f_ext, flm, L, spin, verbosity)

    f(0:L-2, 0:2*L-2) = f_ext(0:L-2, 0:2*L-2)
    f_sp = f_ext(L-1, 0)
    phi_sp = ssht_sampling_mw_p2phi(0, L)

  end subroutine ssht_core_mw_inverse_sp


  !----------------------------------------------------------------------------
  ! ssht_core_mw_forward_sp
  !
  !! South pole wrapper for forward transform using MW method.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [input].
  !!   - f_sp: Function on sphere South pole sample [input].
  !!   - phi_sp: Phi angle corresponding to quoted sample at South pole 
  !!     [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_forward_sp(flm, f, f_sp, phi_sp, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-2 ,0:2*L-2)
    real(dp), intent(in) :: phi_sp
    complex(dpc), intent(in) :: f_sp
    complex(dpc), intent(out) :: flm(0:L**2-1)

    complex(dpc) :: f_ext(0:L-1, 0:2*L-2)
    integer :: p
    real(dp) :: phi
  
    f_ext(0:L-2 ,0:2*L-2) = f(0:L-2 ,0:2*L-2)
    do p = 0, 2*L-2
       phi = ssht_sampling_mw_p2phi(p, L)
       f_ext(L-1, p) = f_sp * exp(I*spin*(phi-phi_sp))
    end do

    call ssht_core_mw_forward_rsp(flm, f_ext, L, spin, verbosity)

  end subroutine ssht_core_mw_forward_sp


  !----------------------------------------------------------------------------
  ! ssht_core_mw_inverse_real_sp
  !
  !! South pole wrapper for inverse transform of real scalar function
  !! using MW method.
  !!
  !! Variables:
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_inverse_real_sp(f, f_sp, flm, L, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: f_sp
    real(dp), intent(out) :: f(0:L-2, 0:2*L-2)

    real(dp) :: f_ext(0:L-1, 0:2*L-2)

    call ssht_core_mw_inverse_real_rsp(f_ext, flm, L, verbosity)

    f(0:L-2, 0:2*L-2) = f_ext(0:L-2, 0:2*L-2)
    f_sp = f_ext(L-1, 0)

  end subroutine ssht_core_mw_inverse_real_sp


  !----------------------------------------------------------------------------
  ! ssht_core_mw_forward_real_sp
  !
  !! South pole wrapper for forward transform of real scalar function
  !! using MW method.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - f(0:L-2, 0:2*L-2): Function on sphere (excluding sample on
  !!     South pole) [output].
  !!   - f_sp: Function on sphere South pole sample [output].
  !!   - L: Harmonic band-limit [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   December 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_forward_real_sp(flm, f, f_sp, L, verbosity)

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    real(dp), intent(in) :: f(0:L-2, 0:2*L-2)
    real(dp), intent(in) :: f_sp
    complex(dpc), intent(out) :: flm(0:L**2-1)

    real(dp) :: f_ext(0:L-1, 0:2*L-2)

    f_ext(0:L-2, 0:2*L-2) = f(0:L-2, 0:2*L-2)
    f_ext(L-1, 0:2*L-2) = f_sp

    call ssht_core_mw_forward_real_rsp(flm, f_ext, L, verbosity)

  end subroutine ssht_core_mw_forward_real_sp


  !============================================================================
  ! Inverse transforms
  !============================================================================


  !----------------------------------------------------------------------------
  ! DH
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_dh_inverse_direct
  !
  !! Compute inverse transform for DH method using direct evaluation
  !! with dlmn evaluated for each theta.
  !!
  !! Notes:
  !!   - Used for code verification purposes (extremely slow!).
  !!
  !! Variables:
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_inverse_direct(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:2*L-1, 0:2*L-2)

    integer :: el, m, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)          
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing inverse transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_inverse_direct...'
       end if
    end if

    f(0:2*L-1 ,0:2*L-2) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do t = 0, 2*L-1
          theta = ssht_sampling_dh_t2theta(t, L)             
          call ssht_dl_beta_operator(dl(-el:el,-el:el), theta, el)
          do m = -el, el
             call ssht_sampling_elm2ind(ind, el, m)
             do p = 0, 2*L-2
                phi = ssht_sampling_dh_p2phi(p, L)
                f(t,p) = f(t,p) + &
                     (-1)**spin * elfactor &
                     * exp(I*m*phi) &
                     * dl(m,-spin) * flm(ind)
             end do
          end do
          
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_inverse_direct


  !----------------------------------------------------------------------------
  ! ssht_core_dh_inverse_direct_factored
  !
  !! Compute inverse transform for DH method using direct evaluation
  !! but dlmn factored, i.e. dlmn evaluated at PI/2.
  !!
  !! Notes:
  !!   - Used for code verification purposes (extremely slow!).
  !!
  !! Variables:
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_inverse_direct_factored(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:2*L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing inverse transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_inverse_direct_factored...'
       end if
    end if     

    f(0:2*L-1 ,0:2*L-2) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             do t = 0, 2*L-1
                theta = ssht_sampling_dh_t2theta(t, L)             
                do p = 0, 2*L-2
                   phi = ssht_sampling_dh_p2phi(p, L)
                   f(t,p) = f(t,p) + &
                        (-1)**spin * elfactor &
                        * exp(-I*PION2*(m+spin)) &
                        * exp(I*m*phi + I*mm*theta) &
                        * dl(mm,m) * dl(mm,-spin) &
                        * flm(ind)
                end do
             end do
          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_inverse_direct_factored


  !----------------------------------------------------------------------------
  ! ssht_core_dh_inverse_sov_direct
  !
  !! Compute inverse transform for DH method using separation of
  !! variables but only discrete (not fast) Fourier transforms.
  !!
  !! Notes:
  !!   - Used for code verification purposes (slow!).
  !!
  !! Variables:
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_inverse_sov_direct(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:2*L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fmt(-(L-1):L-1, 0:2*L-1)

    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing inverse transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_inverse_sov_direct...'
       end if
    end if

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Compute fmt.
    fmt(-(L-1):L-1, 0:2*L-1) = cmplx(0d0, 0d0)
    do t = 0, 2*L-1     
       theta = ssht_sampling_dh_t2theta(t, L)       
       do m = -(L-1), L-1          
          do mm = -(L-1), L-1
             fmt(m,t) = fmt(m,t) + &
                  Fmm(m,mm) * exp(I*mm*theta)
          end do
       end do
    end do

    ! Compute f.
    f(0:2*L-1 ,0:2*L-2) = cmplx(0d0, 0d0)
    do t = 0, 2*L-1       
       do p = 0, 2*L-2
          phi = ssht_sampling_dh_p2phi(p, L)
          do m = -(L-1), L-1          
             f(t,p) = f(t,p) + &
                  fmt(m,t) * exp(I*m*phi)
          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_inverse_sov_direct


  !----------------------------------------------------------------------------
  ! ssht_core_dh_inverse_sov
  !
  !! Compute inverse transform for DH method using separation of
  !! variables and fast Fourier transforms.
  !!
  !! Variables:
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_inverse_sov(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:2*L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fmt(-(L-1):L-1, 0:2*L-1)
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing inverse transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_inverse_sov...'
       end if
    end if


    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Compute fmt.
    fmt(-(L-1):L-1, 0:2*L-1) = cmplx(0d0, 0d0)
    do t = 0, 2*L-1     
       theta = ssht_sampling_dh_t2theta(t, L)       
       do m = -(L-1), L-1          
          do mm = -(L-1), L-1
             fmt(m,t) = fmt(m,t) + &
                  Fmm(m,mm) * exp(I*mm*theta)
          end do
       end do
    end do

    ! Compute f using FFT.
    f(0:2*L-1 ,0:2*L-2) = cmplx(0d0, 0d0)
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, f(0,0:2*L-2), &
         f(0,0:2*L-2), FFTW_BACKWARD, FFTW_MEASURE)
    do t = 0, 2*L-1       

       ! Spatial shift in frequency.
       f(t,0:L-1) = fmt(0:L-1,t)
       f(t,L:2*L-2) = fmt(-(L-1):-1,t)

       call dfftw_execute_dft(fftw_plan, f(t,0:2*L-2), f(t,0:2*L-2))

    end do
    call dfftw_destroy_plan(fftw_plan)

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_inverse_sov


  !----------------------------------------------------------------------------
  ! ssht_core_dh_inverse_sov_sym
  !
  !! Compute inverse transform for DH method using separation of
  !! variables, fast Fourier transforms and exploiting all symmetries
  !! (for complex spin signal).
  !!
  !! Variables:
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_inverse_sov_sym(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:2*L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, 0:L-1)
    complex(dpc) :: fmt(-(L-1):L-1, 0:2*L-1)
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing inverse transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_inverse_sov_sym...'
       end if
    end if

    ! Compute Fmm.
    Fmm(-(L-1):L-1, 0:L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = 0, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Use symmetry to compute Fmm for negative mm.
!!$    do m = -(L-1), L-1       
!!$       do mm = -(L-1), -1
!!$          Fmm(m,mm) = (-1)**(m+spin) * Fmm(m,-mm)
!!$       end do
!!$    end do

    ! Compute fmt.
!!$    fmt(-(L-1):L-1, 0:2*L-1) = cmplx(0d0, 0d0)
!!$    do t = 0, 2*L-1     
!!$       theta = ssht_sampling_dh_t2theta(t, L)       
!!$       do m = -(L-1), L-1          
!!$          do mm = -(L-1), L-1
!!$             fmt(m,t) = fmt(m,t) + &
!!$                  Fmm(m,mm) * exp(I*mm*theta)
!!$          end do
!!$       end do
!!$    end do


    fmt(-(L-1):L-1, 0:2*L-1) = cmplx(0d0, 0d0)
    do t = 0, 2*L-1     
       theta = ssht_sampling_dh_t2theta(t, L)       
       do m = -(L-1), L-1          
          fmt(m,t) = fmt(m,t) + Fmm(m,0)
          do mm = 1, L-1
             fmt(m,t) = fmt(m,t) + &
!!$                  Fmm(m,mm) * exp(I*mm*theta) &
!!$                  + Fmm(m,-mm) * exp(-I*mm*theta)
                  Fmm(m,mm) * (exp(I*mm*theta)  + (-1)**(m+spin) * exp(-I*mm*theta))

          end do
       end do
    end do


    ! Compute f using FFT.
    f(0:2*L-1 ,0:2*L-2) = cmplx(0d0, 0d0)
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, f(0,0:2*L-2), &
         f(0,0:2*L-2), FFTW_BACKWARD, FFTW_MEASURE)
    do t = 0, 2*L-1       

       ! Spatial shift in frequency.
       f(t,0:L-1) = fmt(0:L-1,t)
       f(t,L:2*L-2) = fmt(-(L-1):-1,t)

       call dfftw_execute_dft(fftw_plan, f(t,0:2*L-2), f(t,0:2*L-2))

    end do
    call dfftw_destroy_plan(fftw_plan)

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_inverse_sov_sym


  !----------------------------------------------------------------------------
  ! ssht_core_dh_inverse_sov_sym_real
  !
  !! Compute inverse transform for DH method of real scalar signal
  !! using separation of variables, fast Fourier transforms and
  !! exploiting all symmetries (including additional symmetries for
  !! real signals).
  !!
  !! Variables:
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_inverse_sov_sym_real(f, flm, L, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: f(0:2*L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(0:L-1, 0:L-1)
    complex(dpc) :: Fmm(0:L-1, 0:L-1)
    complex(dpc) :: fmt(0:L-1, 0:2*L-1)
    integer*8 :: fftw_plan

    character(len=STRING_LEN) :: format_spec

    integer :: spin
    integer :: eltmp
    real(dp) :: signs(0:L), ssign
    real(dp) :: dl_mm_spin
    real(dp) :: sqrt_tbl(0:2*(L-1)+1)

    ! Perform precomputations.
    do el = 0, 2*(L-1) + 1
       sqrt_tbl(el) = dsqrt(real(el,kind=dp))
    end do
    do m = 0, L-1, 2
       signs(m)   =  1.0_dp
       signs(m+1) = -1.0_dp
    enddo

    ! Set spin to zero.
    spin = 0
    ssign = signs(spin)    

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing inverse transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',TRUE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_inverse_sov_sym_real...'
       end if
    end if

    ! Compute Fmm.
    Fmm(0:L-1, 0:L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       if (el /= 0 .and. el == abs(spin)) then
          ! Recurse Wigner plane from 0 up to first el, i.e. abs(spin).
          do eltmp = 0, abs(spin)
             call ssht_dl_halfpi_trapani_eighth_table(dl(0:eltmp,0:eltmp), eltmp, &
                  sqrt_tbl(0:2*eltmp+1))
          end do
       else
          call ssht_dl_halfpi_trapani_eighth_table(dl(0:el,0:el), el, &
               sqrt_tbl(0:2*el+1))
          call ssht_dl_halfpi_trapani_fill_eighth2quarter(dl(0:el,0:el), el)
       end if

!       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = 0, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = 0, el
             if (spin <= 0) then
                dl_mm_spin = dl(mm,-spin)
             else
                dl_mm_spin =  signs(el) * signs(mm) * dl(mm,spin)
             end if
             Fmm(m,mm) = Fmm(m,mm) + &
                  ssign * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl_mm_spin &
                  * flm(ind)
          end do
       end do
    end do

    ! Compute fmt.
    fmt(0:L-1, 0:2*L-1) = cmplx(0d0, 0d0)
    do t = 0, 2*L-1     
       theta = ssht_sampling_dh_t2theta(t, L)       
! Switching order of mm and m loops changes computation time by a factor of ~10!
       do m = 0, L-1          
          fmt(m,t) = fmt(m,t) + Fmm(m,0)
       end do
       do mm = 1, L-1
          do m = 0, L-1          
             fmt(m,t) = fmt(m,t) + &
                  Fmm(m,mm) * (exp(I*mm*theta)  + signs(m) * ssign * exp(-I*mm*theta))
          end do
       end do
    end do

    ! Compute f using FFT.
    f(0:2*L-1 ,0:2*L-2) = cmplx(0d0, 0d0)
    call dfftw_plan_dft_c2r_1d(fftw_plan, 2*L-1, Fmm(0:L-1,0), &
         f(0,0:2*L-2), FFTW_MEASURE)
    do t = 0, 2*L-1       
       call dfftw_execute_dft_c2r(fftw_plan, fmt(0:L-1,t), f(t,0:2*L-2))
    end do
    call dfftw_destroy_plan(fftw_plan)

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_inverse_sov_sym_real


  !----------------------------------------------------------------------------
  ! GL
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_gl_inverse_sov_sym
  !
  !! Compute inverse transform for GL method using separation of
  !! variables, fast Fourier transforms and exploiting all symmetries
  !! (for complex spin signal).
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_gl_inverse_sov_sym(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, 0:L-1)
    complex(dpc) :: fmt(-(L-1):L-1, 0:L-1)
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    real(dp) :: thetas(0:L-1)
    real(dp) :: weights(0:L-1)

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing inverse transform using Gauss-Legendre sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_gl_inverse_sov_sym...'
       end if
    end if

    ! Compute weights and theta positions.
    call ssht_sampling_gl_thetas_weights(thetas, weights, L)

    ! Compute Fmm.
    Fmm(-(L-1):L-1, 0:L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = 0, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Use symmetry to compute Fmm for negative mm.
!!$    do m = -(L-1), L-1       
!!$       do mm = -(L-1), -1
!!$          Fmm(m,mm) = (-1)**(m+spin) * Fmm(m,-mm)
!!$       end do
!!$    end do

    ! Compute fmt.
!!$    fmt(-(L-1):L-1, 0:2*L-1) = cmplx(0d0, 0d0)
!!$    do t = 0, 2*L-1     
!!$       theta = ssht_sampling_dh_t2theta(t, L)       
!!$       do m = -(L-1), L-1          
!!$          do mm = -(L-1), L-1
!!$             fmt(m,t) = fmt(m,t) + &
!!$                  Fmm(m,mm) * exp(I*mm*theta)
!!$          end do
!!$       end do
!!$    end do


    fmt(-(L-1):L-1, 0:L-1) = cmplx(0d0, 0d0)
    do t = 0, L-1     
       theta = thetas(t)
       do m = -(L-1), L-1          
          fmt(m,t) = fmt(m,t) + Fmm(m,0)
          do mm = 1, L-1
             fmt(m,t) = fmt(m,t) + &
!!$                  Fmm(m,mm) * exp(I*mm*theta) &
!!$                  + Fmm(m,-mm) * exp(-I*mm*theta)
                  Fmm(m,mm) * (exp(I*mm*theta)  + (-1)**(m+spin) * exp(-I*mm*theta))

          end do
       end do
    end do

    ! Compute f using FFT.
    f(0:L-1 ,0:2*L-2) = cmplx(0d0, 0d0)
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, f(0,0:2*L-2), &
         f(0,0:2*L-2), FFTW_BACKWARD, FFTW_MEASURE)
    do t = 0, L-1       

       ! Spatial shift in frequency.
       f(t,0:L-1) = fmt(0:L-1,t)
       f(t,L:2*L-2) = fmt(-(L-1):-1,t)

       call dfftw_execute_dft(fftw_plan, f(t,0:2*L-2), f(t,0:2*L-2))

    end do
    call dfftw_destroy_plan(fftw_plan)

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_gl_inverse_sov_sym


  !----------------------------------------------------------------------------
  ! ssht_core_gl_inverse_sov_sym_real
  !
  !! Compute inverse transform for GL method of real scalar signal
  !! using separation of variables, fast Fourier transforms and
  !! exploiting all symmetries (including additional symmetries for
  !! real signals).
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_gl_inverse_sov_sym_real(f, flm, L, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(0:L-1, 0:L-1)
    complex(dpc) :: Fmm(0:L-1, 0:L-1)
    complex(dpc) :: fmt(0:L-1, 0:L-1)
    integer*8 :: fftw_plan

    character(len=STRING_LEN) :: format_spec

    real(dp) :: thetas(0:L-1)
    real(dp) :: weights(0:L-1)

    integer :: spin
    integer :: eltmp
    real(dp) :: signs(0:L), ssign
    real(dp) :: dl_mm_spin
    real(dp) :: sqrt_tbl(0:2*(L-1)+1)

    ! Perform precomputations.
    do el = 0, 2*(L-1) + 1
       sqrt_tbl(el) = dsqrt(real(el,kind=dp))
    end do
    do m = 0, L-1, 2
       signs(m)   =  1.0_dp
       signs(m+1) = -1.0_dp
    enddo

    ! Set spin to zero.
    spin = 0
    ssign = signs(spin)
    
    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing inverse transform using Gauss-Legendre sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',TRUE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_gl_inverse_sov_sym_real...'
       end if
    end if

    ! Compute weights and theta positions.
    call ssht_sampling_gl_thetas_weights(thetas, weights, L)

    ! Compute Fmm.
    Fmm(0:L-1, 0:L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
      if (el /= 0 .and. el == abs(spin)) then
          ! Recurse Wigner plane from 0 up to first el, i.e. abs(spin).
          do eltmp = 0, abs(spin)
             call ssht_dl_halfpi_trapani_eighth_table(dl(0:eltmp,0:eltmp), eltmp, &
                  sqrt_tbl(0:2*eltmp+1))
          end do       
       else
          call ssht_dl_halfpi_trapani_eighth_table(dl(0:el,0:el), el, &
               sqrt_tbl(0:2*el+1))
          call ssht_dl_halfpi_trapani_fill_eighth2quarter(dl(0:el,0:el), el)
       end if

       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = 0, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = 0, el
             if (spin <= 0) then
                dl_mm_spin = dl(mm,-spin)
             else
                dl_mm_spin =  signs(el) * signs(mm) * dl(mm,spin)
             end if

             Fmm(m,mm) = Fmm(m,mm) + &
                  ssign * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl_mm_spin &
                  * flm(ind)
          end do
       end do
    end do

    ! Compute fmt.
    fmt(0:L-1, 0:L-1) = cmplx(0d0, 0d0)
    do t = 0, L-1     
       theta = thetas(t)
       do m = 0, L-1          
          fmt(m,t) = fmt(m,t) + Fmm(m,0)
       end do
       do mm = 1, L-1
          do m = 0, L-1          
             fmt(m,t) = fmt(m,t) + &
                  Fmm(m,mm) * (exp(I*mm*theta)  + signs(m) * ssign * exp(-I*mm*theta))
          end do
       end do
    end do

    ! Compute f using FFT.
    f(0:L-1 ,0:2*L-2) = cmplx(0d0, 0d0)
    call dfftw_plan_dft_c2r_1d(fftw_plan, 2*L-1, Fmm(0:L-1,0), &
         f(0,0:2*L-2), FFTW_MEASURE)
    do t = 0, L-1       
       call dfftw_execute_dft_c2r(fftw_plan, fmt(0:L-1,t), f(t,0:2*L-2))
    end do
    call dfftw_destroy_plan(fftw_plan)

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_gl_inverse_sov_sym_real


  !----------------------------------------------------------------------------
  ! MWEO
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_inverse_direct
  !
  !! Compute inverse transform for MWEO method using direct evaluation
  !! with dlmn evaluated for each theta.
  !!
  !! Notes:
  !!   - Used for code verification purposes (extremely slow!).
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_inverse_direct(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing inverse transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_inverse_direct...'
       end if
    end if

    f(0:2*L-1 ,0:2*L-2) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do t = 0, L-1
          theta = ssht_sampling_mweo_t2theta(t, L)             
          call ssht_dl_beta_operator(dl(-el:el,-el:el), theta, el)
          do m = -el, el
             call ssht_sampling_elm2ind(ind, el, m)
             do p = 0, 2*L-2
                phi = ssht_sampling_mweo_p2phi(p, L)
                f(t,p) = f(t,p) + &
                     (-1)**spin * elfactor &
                     * exp(I*m*phi) &
                     * dl(m,-spin) * flm(ind)
             end do
          end do
          
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_inverse_direct


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_inverse_sov_direct
  !
  !! Compute inverse transform for MWEO method using separation of
  !! variables but only discrete (not fast) Fourier transforms.
  !!
  !! Notes:
  !!   - Used for code verification purposes (slow!).
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_inverse_sov_direct(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fext(0:2*L-2, 0:2*L-2)
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing inverse transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_inverse_sov_direct...'
       end if
    end if

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Compute fext using 2D DFT.
    fext(0:2*L-2, 0:2*L-2) = cmplx(0d0, 0d0)
    do t = 0, 2*L-2     
       theta = ssht_sampling_mweo_t2theta(t, L)    
       do p = 0, 2*L-2
          phi = ssht_sampling_mweo_p2phi(p, L)
          do m = -(L-1), L-1          
             do mm = -(L-1), L-1
                fext(t,p) = fext(t,p) + &
                     Fmm(m,mm) * exp(I*m*phi + I*mm*theta)
             end do
          end do
       end do
    end do

    ! Extract f from version of f extended to the torus (fext).
    f(0:L-1, 0:2*L-2) = fext(0:L-1, 0:2*L-2)

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_inverse_sov_direct


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_inverse_sov
  !
  !! Compute inverse transform for MWEO method using separation of
  !! variables and fast Fourier transforms.
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_inverse_sov(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fext(0:2*L-2, 0:2*L-2)
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing inverse transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_inverse_sov...'
       end if
    end if

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Apply phase modulation to account for sampling offset.
!!$    do m = 0, 2*L-2
!!$       do mm = 0, 2*L-2
!!$          Fmm(m-L+1,mm-L+1) = Fmm(m-L+1,mm-L+1) * exp(I*(m+mm)*PI/(2d0*L-1d0))
!!$       end do
!!$    end do
    do m = -(L-1), L-1
       do mm = -(L-1), L-1
          Fmm(m,mm) = Fmm(m,mm) * exp(I*(m+mm)*PI/(2d0*L-1d0))
       end do
    end do  
    
    ! Apply spatial shift.
    fext(0:L-1, 0:L-1) = Fmm(0:L-1,0:L-1)
    fext(L:2*L-2, 0:L-1) = Fmm(-(L-1):-1,0:L-1)
    fext(0:L-1, L:2*L-2) = Fmm(0:L-1,-(L-1):-1)
    fext(L:2*L-2, L:2*L-2) = Fmm(-(L-1):-1,-(L-1):-1)
! If apply phase shift below then just copy Fmm.
!!$    fext(0:2*L-2, 0:2*L-2) = Fmm(-(L-1):L-1,-(L-1):L-1)

    ! Perform 2D FFT.
    ! ** NOTE THAT 2D FFTW SWITCHES DIMENSIONS! HENCE TRANSPOSE BELOW. **
    call dfftw_plan_dft_2d(fftw_plan, 2*L-1, 2*L-1, fext(0:2*L-2,0:2*L-2), &
         fext(0:2*L-2,0:2*L-2), FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(fftw_plan, fext(0:2*L-2,0:2*L-2), fext(0:2*L-2,0:2*L-2))
    call dfftw_destroy_plan(fftw_plan)

    ! Extract f from version of f extended to the torus (fext).
    ! Note transpose due to 2D FFTW (see comment above)!
    ! Also note that additional phase modulation is again applied to 
    ! account for sampling offset.
    f(0:L-1, 0:2*L-2) = transpose(fext(0:2*L-2, 0:L-1)) !&
!         * exp(-I*(L-1)*2d0*PI/(2d0*L-1d0))

! If don't apply spatial shift above then apply phase modulation here.
!!$    do t = 0, L-1     
!!$       theta = ssht_sampling_mweo_t2theta(t, L)    
!!$       do p = 0, 2*L-2
!!$          phi = ssht_sampling_mweo_p2phi(p, L)
!!$          f(t,p) = f(t,p) * exp(-I*(L-1)*(theta+phi))
!!$       end do
!!$    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_inverse_sov


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_inverse_sov_sym
  !
  !! Compute inverse transform for MWEO method using separation of
  !! variables, fast Fourier transforms and exploiting all symmetries
  !! (for complex spin signal).
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_inverse_sov_sym(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fext(0:2*L-2, 0:2*L-2)
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing inverse transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_inverse_sov_sym...'
       end if
    end if

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = 0, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Use symmetry to compute Fmm for negative mm.
    do m = -(L-1), L-1       
       do mm = -(L-1), -1
          Fmm(m,mm) = (-1)**(m+spin) * Fmm(m,-mm)
       end do
    end do

    ! Apply phase modulation to account for sampling offset.
!!$    do m = 0, 2*L-2
!!$       do mm = 0, 2*L-2
!!$          Fmm(m-L+1,mm-L+1) = Fmm(m-L+1,mm-L+1) * exp(I*(m+mm)*PI/(2d0*L-1d0))
!!$       end do
!!$    end do
    do m = -(L-1), L-1
       do mm = -(L-1), L-1
          Fmm(m,mm) = Fmm(m,mm) * exp(I*(m+mm)*PI/(2d0*L-1d0))
       end do
    end do    


    ! Apply spatial shift.
    fext(0:L-1, 0:L-1) = Fmm(0:L-1,0:L-1)
    fext(L:2*L-2, 0:L-1) = Fmm(-(L-1):-1,0:L-1)
    fext(0:L-1, L:2*L-2) = Fmm(0:L-1,-(L-1):-1)
    fext(L:2*L-2, L:2*L-2) = Fmm(-(L-1):-1,-(L-1):-1)
! If apply phase shift below then just copy Fmm.
!!$    fext(0:2*L-2, 0:2*L-2) = Fmm(-(L-1):L-1,-(L-1):L-1)

    ! Perform 2D FFT.
    ! ** NOTE THAT 2D FFTW SWITCHES DIMENSIONS! HENCE TRANSPOSE BELOW. **
    call dfftw_plan_dft_2d(fftw_plan, 2*L-1, 2*L-1, fext(0:2*L-2,0:2*L-2), &
         fext(0:2*L-2,0:2*L-2), FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(fftw_plan, fext(0:2*L-2,0:2*L-2), fext(0:2*L-2,0:2*L-2))
    call dfftw_destroy_plan(fftw_plan)

    ! Extract f from version of f extended to the torus (fext).
    ! Note transpose due to 2D FFTW (see comment above)!
    ! Also note that additional phase modulation is again applied to 
    ! account for sampling offset.
    f(0:L-1, 0:2*L-2) = transpose(fext(0:2*L-2, 0:L-1))! &
!         * exp(-I*(L-1)*2d0*PI/(2d0*L-1d0))

! If don't apply spatial shift above then apply phase modulation here.
!!$    do t = 0, L-1     
!!$       theta = ssht_sampling_mweo_t2theta(t, L)    
!!$       do p = 0, 2*L-2
!!$          phi = ssht_sampling_mweo_p2phi(p, L)
!!$          f(t,p) = f(t,p) * exp(-I*(L-1)*(theta+phi))
!!$       end do
!!$    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_inverse_sov_sym


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_inverse_sov_sym_real
  !
  !! Compute inverse transform for MWEO method of real scalar signal
  !! using separation of variables, fast Fourier transforms and
  !! exploiting all symmetries (including additional symmetries for
  !! real signals).
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_inverse_sov_sym_real(f, flm, L, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(0:L-1, -(L-1):L-1)
    complex(dpc) :: fext(0:2*L-2, 0:2*L-2)
    real(dp) :: fext_real(0:2*L-2,0:2*L-2)
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    integer :: spin

    spin = 0

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing inverse transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',TRUE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_inverse_sov_sym_real...'
       end if
    end if

    ! Compute Fmm.
    Fmm(0:L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = 0, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = 0, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Use symmetry to compute Fmm for negative mm.
    do m = 0, L-1       
       do mm = -(L-1), -1
          Fmm(m,mm) = (-1)**(m+spin) * Fmm(m,-mm)
       end do
    end do

    ! Apply phase modulation to account for sampling offset.
    do m = 0, L-1
       do mm = -(L-1), L-1
          Fmm(m,mm) = Fmm(m,mm) * exp(I*(m+mm)*PI/(2d0*L-1d0))
       end do
    end do
    
    ! Apply spatial shift.
    fext(0:L-1, 0:L-1) = Fmm(0:L-1,0:L-1)
    fext(0:L-1, L:2*L-2) = Fmm(0:L-1,-(L-1):-1)

    ! Perform 2D FFT.
    ! ** NOTE THAT 2D FFTW SWITCHES DIMENSIONS! HENCE TRANSPOSE BELOW. **
    call dfftw_plan_dft_c2r_2d(fftw_plan, 2*L-1, 2*L-1, fext(0:L-1,0:2*L-2), &
         fext_real(0:2*L-2,0:2*L-2), FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(fftw_plan, fext(0:L-1,0:2*L-2), fext_real(0:2*L-2,0:2*L-2))
    call dfftw_destroy_plan(fftw_plan)

    ! Extract f from version of f extended to the torus (fext).
    ! Note transpose due to 2D FFTW (see comment above)!
    ! Also note that additional phase modulation is again applied to 
    ! account for sampling offset.
    f(0:L-1, 0:2*L-2) = transpose(fext_real(0:2*L-2, 0:L-1))

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_inverse_sov_sym_real


  !----------------------------------------------------------------------------
  ! MW
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_mw_inverse_sov_direct
  !
  !! Compute inverse transform for MW method using separation of
  !! variables but only discrete (not fast) Fourier transforms.
  !!
  !! Notes:
  !!   - Used for code verification purposes (slow!).
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_inverse_sov_direct(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fext(0:2*L-2, 0:2*L-2)
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing inverse transform using McEwen and Wiaux', &
               ' sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_inverse_direct...'
       end if
    end if

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Compute fext using 2D DFT.
    fext(0:2*L-2, 0:2*L-2) = cmplx(0d0, 0d0)
    do t = 0, 2*L-2     
       theta = ssht_sampling_mw_t2theta(t, L)    
       do p = 0, 2*L-2
          phi = ssht_sampling_mw_p2phi(p, L)
          do m = -(L-1), L-1          
             do mm = -(L-1), L-1
                fext(t,p) = fext(t,p) + &
                     Fmm(m,mm) * exp(I*m*phi + I*mm*theta)
             end do
          end do
       end do
    end do

    ! Extract f from version of f extended to the torus (fext).
    f(0:L-1, 0:2*L-2) = fext(0:L-1, 0:2*L-2)

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_mw_inverse_sov_direct


  !----------------------------------------------------------------------------
  ! ssht_core_mw_inverse_sov
  !
  !! Compute inverse transform for MW method using separation of
  !! variables and fast Fourier transforms.
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_inverse_sov(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fext(0:2*L-2, 0:2*L-2)
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing inverse transform using McEwen and Wiaux', &
               ' sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mw_inverse_sov...'
       end if
    end if

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Apply phase modulation to account for sampling offset.
    do mm = -(L-1), L-1
       Fmm(-(L-1):L-1,mm) = Fmm(-(L-1):L-1,mm) * exp(I*mm*PI/(2d0*L-1d0))
    end do

    ! Apply spatial shift.
    fext(0:L-1, 0:L-1) = Fmm(0:L-1,0:L-1)
    fext(L:2*L-2, 0:L-1) = Fmm(-(L-1):-1,0:L-1)
    fext(0:L-1, L:2*L-2) = Fmm(0:L-1,-(L-1):-1)
    fext(L:2*L-2, L:2*L-2) = Fmm(-(L-1):-1,-(L-1):-1)

    ! Perform 2D FFT.
    ! ** NOTE THAT 2D FFTW SWITCHES DIMENSIONS! HENCE TRANSPOSE BELOW. **
    call dfftw_plan_dft_2d(fftw_plan, 2*L-1, 2*L-1, fext(0:2*L-2,0:2*L-2), &
         fext(0:2*L-2,0:2*L-2), FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(fftw_plan, fext(0:2*L-2,0:2*L-2), fext(0:2*L-2,0:2*L-2))
    call dfftw_destroy_plan(fftw_plan)

    ! Extract f from version of f extended to the torus (fext).
    f(0:L-1, 0:2*L-2) = transpose(fext(0:2*L-2, 0:L-1))

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_mw_inverse_sov


  !----------------------------------------------------------------------------
  ! ssht_core_mw_inverse_sov_sym
  !
  !! Compute inverse transform for MW method using separation of
  !! variables, fast Fourier transforms and exploiting all symmetries
  !! (for complex spin signal).
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_inverse_sov_sym(f, flm, L, spin, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    complex(dpc), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fext(0:2*L-2, 0:2*L-2)
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing inverse transform using McEwen and Wiaux', &
               ' sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mw_inverse_sov_sym...'
       end if
    end if

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = 0, el
             Fmm(m,mm) = Fmm(m,mm) + &
                  (-1)**spin * elfactor &
                  * exp(-I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * flm(ind)
          end do
       end do
    end do

    ! Use symmetry to compute Fmm for negative mm.
    do m = -(L-1), L-1       
       do mm = -(L-1), -1
          Fmm(m,mm) = (-1)**(m+spin) * Fmm(m,-mm)
       end do
    end do

    ! Apply phase modulation to account for sampling offset.
    do mm = -(L-1), L-1
       Fmm(-(L-1):L-1,mm) = Fmm(-(L-1):L-1,mm) * exp(I*mm*PI/(2d0*L-1d0))
    end do

    ! Apply spatial shift.
    fext(0:L-1, 0:L-1) = Fmm(0:L-1,0:L-1)
    fext(L:2*L-2, 0:L-1) = Fmm(-(L-1):-1,0:L-1)
    fext(0:L-1, L:2*L-2) = Fmm(0:L-1,-(L-1):-1)
    fext(L:2*L-2, L:2*L-2) = Fmm(-(L-1):-1,-(L-1):-1)

    ! Perform 2D FFT.
    ! ** NOTE THAT 2D FFTW SWITCHES DIMENSIONS! HENCE TRANSPOSE BELOW. **
    call dfftw_plan_dft_2d(fftw_plan, 2*L-1, 2*L-1, fext(0:2*L-2,0:2*L-2), &
         fext(0:2*L-2,0:2*L-2), FFTW_BACKWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(fftw_plan, fext(0:2*L-2,0:2*L-2), fext(0:2*L-2,0:2*L-2))
    call dfftw_destroy_plan(fftw_plan)

    ! Extract f from version of f extended to the torus (fext).
    f(0:L-1, 0:2*L-2) = transpose(fext(0:2*L-2, 0:L-1))

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_mw_inverse_sov_sym


  !----------------------------------------------------------------------------
  ! ssht_core_mw_inverse_sov_sym_real
  !
  !! Compute inverse transform for MW method of real scalar signal
  !! using separation of variables, fast Fourier transforms and
  !! exploiting all symmetries (including additional symmetries for
  !! real signals).
  !!
  !! Variables:
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [output].
  !!   - flm(0:L**2-1): Harmonic coefficients [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_inverse_sov_sym_real(f, flm, L, verbosity)
    
    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: flm(0:L**2-1)
    real(dp), intent(out) :: f(0:L-1, 0:2*L-2)

    integer :: el, m, mm, t, p, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: dl(0:L-1, 0:L-1)
    complex(dpc) :: Fmm(0:L-1, -(L-1):L-1)
    complex(dpc) :: fext(0:L-1, 0:2*L-2)
    real(dp) :: fext_real(0:2*L-2,0:2*L-2)
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    integer :: spin
    integer :: eltmp
    real(dp) :: signs(0:L), ssign
    real(dp) :: dl_mm_spin
    real(dp) :: sqrt_tbl(0:2*(L-1)+1)

    ! Perform precomputations.
    do el = 0, 2*(L-1) + 1
       sqrt_tbl(el) = dsqrt(real(el,kind=dp))
    end do
    do m = 0, L-1, 2
       signs(m)   =  1.0_dp
       signs(m+1) = -1.0_dp
    enddo

    ! Set spin to zero.
    spin = 0
    ssign = signs(spin)

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing inverse transform using McEwen and Wiaux', &
               ' sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',TRUE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mw_inverse_sov_sym...'
       end if
    end if

    ! Compute Fmm.
    Fmm(0:L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       if (el /= 0 .and. el == abs(spin)) then
          ! Recurse Wigner plane from 0 up to first el, i.e. abs(spin).
          do eltmp = 0, abs(spin)
             call ssht_dl_halfpi_trapani_eighth_table(dl(0:eltmp,0:eltmp), eltmp, &
                  sqrt_tbl(0:2*eltmp+1))
          end do       
       else
          call ssht_dl_halfpi_trapani_eighth_table(dl(0:el,0:el), el, &
               sqrt_tbl(0:2*el+1))
          call ssht_dl_halfpi_trapani_fill_eighth2quarter(dl(0:el,0:el), el)
       end if

       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = 0, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = 0, el
             if (spin <= 0) then
                dl_mm_spin = dl(mm,-spin)
             else
                dl_mm_spin =  signs(el) * signs(mm) * dl(mm,spin)
             end if

             Fmm(m,mm) = Fmm(m,mm) + &
                  ssign * elfactor &
                  * exp(-I*PION2*(m+spin)) &
!!$                  * dl(mm,m) * dl(mm,-spin) &
                  * dl(mm,m) * dl_mm_spin &
                  * flm(ind)
          end do
       end do
    end do

    ! Use symmetry to compute Fmm for negative mm.
    do m = 0, L-1       
       do mm = -(L-1), -1
          Fmm(m,mm) = signs(m) * ssign * Fmm(m,-mm)
       end do
    end do

    ! Apply phase modulation to account for sampling offset.
    do mm = -(L-1), L-1
       Fmm(0:L-1,mm) = Fmm(0:L-1,mm) * exp(I*mm*PI/(2d0*L-1d0))
    end do

    ! Apply spatial shift.
    fext(0:L-1, 0:L-1) = Fmm(0:L-1,0:L-1)
    fext(0:L-1, L:2*L-2) = Fmm(0:L-1,-(L-1):-1)

    ! Perform 2D FFT.
    ! ** NOTE THAT 2D FFTW SWITCHES DIMENSIONS! HENCE TRANSPOSE BELOW. **
    call dfftw_plan_dft_c2r_2d(fftw_plan, 2*L-1, 2*L-1, fext(0:L-1,0:2*L-2), &
         fext_real(0:2*L-2,0:2*L-2), FFTW_ESTIMATE)
    call dfftw_execute_dft_c2r(fftw_plan, fext(0:L-1,0:2*L-2), fext_real(0:2*L-2,0:2*L-2))
    call dfftw_destroy_plan(fftw_plan)

    ! Extract f from version of f extended to the torus (fext).
    f(0:L-1, 0:2*L-2) = transpose(fext_real(0:2*L-2, 0:L-1))

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Inverse transform computed!'
       end if
    end if

  end subroutine ssht_core_mw_inverse_sov_sym_real


  !============================================================================
  ! Forward transforms
  !============================================================================


  !----------------------------------------------------------------------------
  ! DH
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_dh_forward_sov_direct
  !
  !! Compute forward transform for DH method using separation of
  !! variables but only discrete (not fast) Fourier transforms.
  !!
  !! Notes:
  !!   - Used for code verification purposes (extremely slow!).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_forward_sov_direct(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:2*L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: w
    complex(dpc) :: fmt(-(L-1):L-1, 0:2*L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing forward transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_forward_sov_direct...'
       end if
    end if

    ! Compute fmt.
    fmt(-(L-1):L-1, 0:2*L-1) = cmplx(0d0, 0d0)
    do p = 0, 2*L-2
       phi = ssht_sampling_dh_p2phi(p, L)
       do m = -(L-1), L-1
          do t = 0, 2*L-1             
             fmt(m,t) = fmt(m,t) + &
                  f(t,p) * exp(-I*m*phi)
          end do
       end do
    end do
    fmt(-(L-1):L-1, 0:2*L-1) = fmt(-(L-1):L-1, 0:2*L-1) &
         * 2d0*PI / (2d0*L-1d0)
    
    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do t = 0, 2*L-1
       theta = ssht_sampling_dh_t2theta(t, L)
       w = ssht_sampling_weight_dh(theta, L)
       do m = -(L-1), L-1
          do mm = -(L-1), L-1
             Fmm(m,mm) = Fmm(m,mm) + &
                  fmt(m,t) * exp(-I*mm*theta) * w
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Fmm(m,mm)
          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_forward_sov_direct


  !----------------------------------------------------------------------------
  ! ssht_core_dh_forward_sov
  !
  !! Compute forward transform for DH method using separation of
  !! variables and fast Fourier transforms.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_forward_sov(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:2*L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: w
    complex(dpc) :: fmt(-(L-1):L-1, 0:2*L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    integer*8 :: fftw_plan
    complex(dpc) :: tmp(0:2*L-2)
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing forward transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_forward_sov...'
       end if
    end if

    ! Compute fmt using FFT.     
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, fmt(-(L-1):L-1,0), &
         fmt(-(L-1):L-1,0), FFTW_FORWARD, FFTW_MEASURE)
    do t = 0, 2*L-1             

       call dfftw_execute_dft(fftw_plan, f(t,0:2*L-2), tmp(0:2*L-2))

       ! Spatial shift in frequency.
       fmt(0:L-1, t) = tmp(0:L-1)
       fmt(-(L-1):-1, t) = tmp(L:2*L-2)
    end do
    call dfftw_destroy_plan(fftw_plan)
    fmt(-(L-1):L-1, 0:2*L-1) = fmt(-(L-1):L-1, 0:2*L-1) &
         * 2d0*PI / (2d0*L-1d0)

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do t = 0, 2*L-1
       theta = ssht_sampling_dh_t2theta(t, L)
       w = ssht_sampling_weight_dh(theta, L)
       do m = -(L-1), L-1
          do mm = -(L-1), L-1
             Fmm(m,mm) = Fmm(m,mm) + &
                  fmt(m,t) * exp(-I*mm*theta) * w
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Fmm(m,mm)
          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_forward_sov


  !----------------------------------------------------------------------------
  ! ssht_core_dh_forward_sov_sym
  !
  !! Compute forward transform for DH method using separation of
  !! variables, fast Fourier transforms and exploiting all symmetries
  !! (for complex spin signal).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_forward_sov_sym(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:2*L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: w
    complex(dpc) :: fmt(-(L-1):L-1, 0:2*L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    integer*8 :: fftw_plan
    complex(dpc) :: tmp(0:2*L-2)
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing forward transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_forward_sov_sym...'
       end if
    end if

    ! Compute fmt using FFT.     
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, fmt(-(L-1):L-1,0), &
         fmt(-(L-1):L-1,0), FFTW_FORWARD, FFTW_MEASURE)
    do t = 0, 2*L-1             

       call dfftw_execute_dft(fftw_plan, f(t,0:2*L-2), tmp(0:2*L-2))

       ! Spatial shift in frequency.
       fmt(0:L-1, t) = tmp(0:L-1)
       fmt(-(L-1):-1, t) = tmp(L:2*L-2)
    end do
    call dfftw_destroy_plan(fftw_plan)
    fmt(-(L-1):L-1, 0:2*L-1) = fmt(-(L-1):L-1, 0:2*L-1) &
         * 2d0*PI / (2d0*L-1d0)

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do t = 0, 2*L-1
       theta = ssht_sampling_dh_t2theta(t, L)
       w = ssht_sampling_weight_dh(theta, L)
       do m = -(L-1), L-1
          do mm = -(L-1), L-1
             Fmm(m,mm) = Fmm(m,mm) + &
                  fmt(m,t) * exp(-I*mm*theta) * w
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)

          flm(ind) = flm(ind) + &
               (-1)**spin * elfactor &
               * exp(I*PION2*(m+spin)) &
               * dl(0,m) * dl(0,-spin) &
               * Fmm(m,0)

          do mm = 1, el
             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * (Fmm(m,mm) + (-1)**(m+spin)*Fmm(m,-mm))

          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_forward_sov_sym


  !----------------------------------------------------------------------------
  ! ssht_core_dh_forward_sov_sym_real
  !
  !! Compute forward transform for DH method of real scalar signal
  !! using separation of variables, fast Fourier transforms and
  !! exploiting all symmetries (including additional symmetries for
  !! real signals).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:2*L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_dh_forward_sov_sym_real(flm, f, L, verbosity)

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    real(dp), intent(in) :: f(0:2*L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: w
    complex(dpc) :: fmt(0:L-1, 0:2*L-1)
    complex(dpc) :: Fmm(0:L-1, -(L-1):L-1)
    real(dp) :: dl(0:L-1, 0:L-1)
    integer*8 :: fftw_plan

    integer :: spin
    real(dp) :: tmp(0:2*L-2)
    integer :: ind_nm
    character(len=STRING_LEN) :: format_spec
    integer :: eltmp

    real(dp) :: signs(0:L), ssign
    real(dp) :: dl_mm_spin, dl_0_spin
    real(dp) :: sqrt_tbl(0:2*(L-1)+1)

    ! Perform precomputations.
    do el = 0, 2*(L-1) + 1
       sqrt_tbl(el) = dsqrt(real(el,kind=dp))
    end do
    do m = 0, L-1, 2
       signs(m)   =  1.0_dp
       signs(m+1) = -1.0_dp
    enddo

    ! Set spin to zero.
    spin = 0
    ssign = signs(spin)

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing forward transform using Driscoll and Healy sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',TRUE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_dh_forward_sov_sym_real...'
       end if
    end if

    ! Compute fmt using FFT.     
    call dfftw_plan_dft_r2c_1d(fftw_plan, 2*L-1, tmp(0:2*L-2), &
         fmt(0:L-1,0), FFTW_MEASURE)
    do t = 0, 2*L-1             
       call dfftw_execute_dft_r2c(fftw_plan, f(t,0:2*L-2), fmt(0:L-1,t))
    end do
    call dfftw_destroy_plan(fftw_plan)
    fmt(0:L-1, 0:2*L-1) = fmt(0:L-1, 0:2*L-1) &
         * 2d0*PI / (2d0*L-1d0)

    ! Compute Fmm.
    Fmm(0:L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do t = 0, 2*L-1
       theta = ssht_sampling_dh_t2theta(t, L)
       w = ssht_sampling_weight_dh(theta, L)
! Switching order of mm and m loops changes computation time by a factor of ~10!
       do mm = -(L-1), L-1
          do m = 0, L-1
             Fmm(m,mm) = Fmm(m,mm) + &
                  fmt(m,t) * exp(-I*mm*theta) * w
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       if (el /= 0 .and. el == abs(spin)) then
          ! Recurse Wigner plane from 0 up to first el, i.e. abs(spin).
          do eltmp = 0, abs(spin)
             call ssht_dl_halfpi_trapani_eighth_table(dl(0:eltmp,0:eltmp), eltmp, &
                  sqrt_tbl(0:2*eltmp+1))
          end do
       else
          call ssht_dl_halfpi_trapani_eighth_table(dl(0:el,0:el), el, &
               sqrt_tbl(0:2*el+1))
          call ssht_dl_halfpi_trapani_fill_eighth2quarter(dl(0:el,0:el), el)
       end if
      
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = 0, el
          call ssht_sampling_elm2ind(ind, el, m)

          if (spin <= 0) then
             dl_0_spin = dl(0,-spin)
          else
             dl_0_spin = signs(el) * dl(0,spin)
          end if

          flm(ind) = flm(ind) + &
               (-1)**spin * elfactor &
               * exp(I*PION2*(m+spin)) &
               * dl(0,m) * dl_0_spin &
               * Fmm(m,0)

          do mm = 1, el

             if (spin <= 0) then
                dl_mm_spin = dl(mm,-spin)
             else
                dl_mm_spin = signs(el) * signs(mm) * dl(mm,spin)
             end if

             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl_mm_spin &
                  * (Fmm(m,mm) + signs(m) * ssign * Fmm(m,-mm))

          end do
       end do
    end do

    ! Set flm values for negative m using conjugate symmetry.
    do el = abs(spin), L-1
       do m = 1, el
          call ssht_sampling_elm2ind(ind, el, m)
          call ssht_sampling_elm2ind(ind_nm, el, -m)
          flm(ind_nm) = (-1)**m * conjg(flm(ind))
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_dh_forward_sov_sym_real


  !----------------------------------------------------------------------------
  ! GL
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_gl_forward_sov_sym
  !
  !! Compute forward transform for GL method using separation of
  !! variables, fast Fourier transforms and exploiting all symmetries
  !! (for complex spin signal).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_gl_forward_sov_sym(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: w
    complex(dpc) :: fmt(-(L-1):L-1, 0:L-1)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    integer*8 :: fftw_plan
    complex(dpc) :: tmp(0:2*L-2)
    character(len=STRING_LEN) :: format_spec

    real(dp) :: thetas(0:L-1)
    real(dp) :: weights(0:L-1)

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing forward transform using Gauss-Legendre sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_gl_forward_sov_sym...'
       end if
    end if

    ! Compute weights and theta positions.
    call ssht_sampling_gl_thetas_weights(thetas, weights, L)

    ! Compute fmt using FFT.     
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, fmt(-(L-1):L-1,0), &
         fmt(-(L-1):L-1,0), FFTW_FORWARD, FFTW_MEASURE)
    do t = 0, L-1             

       call dfftw_execute_dft(fftw_plan, f(t,0:2*L-2), tmp(0:2*L-2))

       ! Spatial shift in frequency.
       fmt(0:L-1, t) = tmp(0:L-1)
       fmt(-(L-1):-1, t) = tmp(L:2*L-2)
    end do
    call dfftw_destroy_plan(fftw_plan)
    fmt(-(L-1):L-1, 0:L-1) = fmt(-(L-1):L-1, 0:L-1) &
         * 2d0*PI / (2d0*L-1d0)

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do t = 0, L-1
       theta = thetas(t)
       w = weights(t)
       do m = -(L-1), L-1
          do mm = -(L-1), L-1
             Fmm(m,mm) = Fmm(m,mm) + &
                  fmt(m,t) * exp(-I*mm*theta) * w
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)

          flm(ind) = flm(ind) + &
               (-1)**spin * elfactor &
               * exp(I*PION2*(m+spin)) &
               * dl(0,m) * dl(0,-spin) &
               * Fmm(m,0)

          do mm = 1, el
             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * (Fmm(m,mm) + (-1)**(m+spin)*Fmm(m,-mm))

          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_gl_forward_sov_sym


  !----------------------------------------------------------------------------
  ! ssht_core_gl_forward_sov_sym_real
  !
  !! Compute forward transform for GL method of real scalar signal
  !! using separation of variables, fast Fourier transforms and
  !! exploiting all symmetries (including additional symmetries for
  !! real signals).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_gl_forward_sov_sym_real(flm, f, L, verbosity)

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    real(dp), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind
    real(dp) :: theta, phi
    real(dp) :: elfactor
    real(dp) :: w
    complex(dpc) :: fmt(0:L-1, 0:L-1)
    complex(dpc) :: Fmm(0:L-1, -(L-1):L-1)
    real(dp) :: dl(0:L-1, 0:L-1)
    integer*8 :: fftw_plan

    integer :: spin
    real(dp) :: tmp(0:2*L-2)
    integer :: ind_nm
    character(len=STRING_LEN) :: format_spec

    real(dp) :: thetas(0:L-1)
    real(dp) :: weights(0:L-1)

    real(dp) :: signs(0:L), ssign
    real(dp) :: dl_mm_spin, dl_0_spin
    real(dp) :: sqrt_tbl(0:2*(L-1)+1)
    integer :: eltmp

    ! Perform precomputations.
    do el = 0, 2*(L-1) + 1
       sqrt_tbl(el) = dsqrt(real(el,kind=dp))
    end do
    do m = 0, L-1, 2
       signs(m)   =  1.0_dp
       signs(m+1) = -1.0_dp
    enddo

    ! Set spin to zero.
    spin = 0
    ssign = signs(spin)

    spin = 0

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Computing forward transform using Gauss-Legendre sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',TRUE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_gl_forward_sov_sym_real...'
       end if
    end if

    ! Compute weights and theta positions.
    call ssht_sampling_gl_thetas_weights(thetas, weights, L)

    ! Compute fmt using FFT.     
    call dfftw_plan_dft_r2c_1d(fftw_plan, 2*L-1, tmp(0:2*L-2), &
         fmt(0:L-1,0), FFTW_MEASURE)
    do t = 0, L-1             
       call dfftw_execute_dft_r2c(fftw_plan, f(t,0:2*L-2), fmt(0:L-1,t))
    end do
    call dfftw_destroy_plan(fftw_plan)
    fmt(0:L-1, 0:L-1) = fmt(0:L-1, 0:L-1) &
         * 2d0*PI / (2d0*L-1d0)

    ! Compute Fmm.
    Fmm(0:L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do t = 0, L-1
       theta = thetas(t)
       w = weights(t)
       do mm = -(L-1), L-1
          do m = 0, L-1
             Fmm(m,mm) = Fmm(m,mm) + &
                  fmt(m,t) * exp(-I*mm*theta) * w
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       if (el /= 0 .and. el == abs(spin)) then
          ! Recurse Wigner plane from 0 up to first el, i.e. abs(spin).
          do eltmp = 0, abs(spin)
             call ssht_dl_halfpi_trapani_eighth_table(dl(0:eltmp,0:eltmp), eltmp, &
                  sqrt_tbl(0:2*eltmp+1))
          end do
       else
          call ssht_dl_halfpi_trapani_eighth_table(dl(0:el,0:el), el, &
               sqrt_tbl(0:2*el+1))
          call ssht_dl_halfpi_trapani_fill_eighth2quarter(dl(0:el,0:el), el)
       end if

       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = 0, el
          call ssht_sampling_elm2ind(ind, el, m)

          if (spin <= 0) then
             dl_0_spin = dl(0,-spin)
          else
             dl_0_spin = signs(el) * dl(0,spin)
          end if         

          flm(ind) = flm(ind) + &
               ssign * elfactor &
               * exp(I*PION2*(m+spin)) &
               * dl(0,m) * dl_0_spin &
               * Fmm(m,0)

          do mm = 1, el

             if (spin <= 0) then
                dl_mm_spin = dl(mm,-spin)
             else
                dl_mm_spin = signs(el) * signs(mm) * dl(mm,spin)
             end if

             flm(ind) = flm(ind) + &
                  ssign * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl_mm_spin &
                  * (Fmm(m,mm) + signs(m) * ssign * Fmm(m,-mm))

          end do
       end do
    end do

    ! Set flm values for negative m using conjugate symmetry.
    do el = abs(spin), L-1
       do m = 1, el
          call ssht_sampling_elm2ind(ind, el, m)
          call ssht_sampling_elm2ind(ind_nm, el, -m)
          flm(ind_nm) = signs(m) * conjg(flm(ind))
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_gl_forward_sov_sym_real


  !----------------------------------------------------------------------------
  ! MWEO
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_forward_sov_direct
  !
  !! Compute forward transform for MWEO method using separation of
  !! variables but only discrete (not fast) Fourier transforms.
  !!
  !! Notes:
  !!   - Used for code verification purposes (extremely slow!).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_forward_sov_direct(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fe(0:2*L-2 ,0:2*L-2)
    complex(dpc) :: fo(0:2*L-2 ,0:2*L-2)
    complex(dpc) :: Fmme(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmmo(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmme(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmmo(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmm_term
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_forward_sov_direct...'
       end if
    end if

    ! Extend f to the torus with even and odd extensions 
    ! about theta=PI.
    fe(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fe(L:2*L-2, 0:2*L-2) = f(L-2:0:-1, 0:2*L-2)
    fo(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fo(L:2*L-2, 0:2*L-2) = -f(L-2:0:-1, 0:2*L-2)

    ! Compute Fmm for even and odd extensions.
    Fmme(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    Fmmo(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do p = 0, 2*L-2
       phi = ssht_sampling_mweo_p2phi(p, L)
       do t = 0, 2*L-2
          theta = ssht_sampling_mw_t2theta(t, L)    
          do m = -(L-1), L-1
             do mm = -(L-1), L-1
                Fmme(m,mm) = Fmme(m,mm) + &
                     fe(t,p) * exp(-I*(m*phi + mm*theta))
                Fmmo(m,mm) = Fmmo(m,mm) + &
                     fo(t,p) * exp(-I*(m*phi + mm*theta))
             end do
          end do

       end do
    end do
    Fmme(-(L-1):L-1, -(L-1):L-1) = Fmme(-(L-1):L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2
    Fmmo(-(L-1):L-1, -(L-1):L-1) = Fmmo(-(L-1):L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2

    ! Compute Gmm for even and odd extensions by direct calculation 
    ! of convolution.
    Gmme(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    Gmmo(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do m = -(L-1), L-1
       do mm = -(L-1), L-1
          do k = -(L-1), L-1 
             Gmme(m,mm) = Gmme(m,mm) + &
                  Fmme(m,k) * ssht_sampling_weight_mw(k - mm) * 2d0 * PI
             Gmmo(m,mm) = Gmmo(m,mm) + &
                  Fmmo(m,k) *ssht_sampling_weight_mw(k - mm) * 2d0 * PI
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             if (mod(m+spin,2) == 0) then
                ! m+spin even
                Gmm_term = Gmme(m,mm)
             else
                ! m+spin odd
                Gmm_term = Gmmo(m,mm)
             end if

             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Gmm_term

          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_forward_sov_direct


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_forward_sov
  !
  !! Compute forward transform for MWEO method using separation of
  !! variables and fast Fourier transforms.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_forward_sov(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fe(0:2*L-2 ,0:2*L-2)
    complex(dpc) :: fo(0:2*L-2 ,0:2*L-2)
    complex(dpc) :: Fmme(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmmo(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmme(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmmo(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmm_term
    integer*8 :: fftw_plan
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_forward_sov...'
       end if
    end if

    ! Extend f to the torus with even and odd extensions 
    ! about theta=PI.
    fe(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fe(L:2*L-2, 0:2*L-2) = f(L-2:0:-1, 0:2*L-2)
    fo(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fo(L:2*L-2, 0:2*L-2) = -f(L-2:0:-1, 0:2*L-2)

! If don't apply spatial shift below then apply phase modulation here.
!!$    do p = 0, 2*L-2
!!$       phi = ssht_sampling_mweo_p2phi(p, L)
!!$       do t = 0, 2*L-2
!!$          theta = ssht_sampling_mweo_t2theta(t, L)   
!!$          fe(t,p) = fe(t,p) * exp(I*(phi+theta)*(L-1))
!!$          fo(t,p) = fo(t,p) * exp(I*(phi+theta)*(L-1))
!!$       end do
!!$    end do

    ! Compute Fmm for even and odd extensions by 2D FFT.
    ! ** NOTE THAT 2D FFTW SWITCHES DIMENSIONS! HENCE TRANSPOSE BELOW. **
    call dfftw_plan_dft_2d(fftw_plan, 2*L-1, 2*L-1, fe(0:2*L-2,0:2*L-2), &
         fe(0:2*L-2,0:2*L-2), FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(fftw_plan, fe(0:2*L-2,0:2*L-2), fe(0:2*L-2,0:2*L-2))
    call dfftw_execute_dft(fftw_plan, fo(0:2*L-2,0:2*L-2), fo(0:2*L-2,0:2*L-2))
    call dfftw_destroy_plan(fftw_plan)

! If apply phase shift above just copy Fmm (and transpose 
! to account for 2D FFT).
!!$    Fmme(-(L-1):L-1, -(L-1):L-1) = transpose(fe(0:2*L-2,0:2*L-2))
!!$    Fmmo(-(L-1):L-1, -(L-1):L-1) = transpose(fo(0:2*L-2,0:2*L-2))

    ! Apply spatial shift in frequency.
    fe(0:2*L-2,0:2*L-2) = transpose(fe(0:2*L-2,0:2*L-2)) * exp(I*(L-1)*2d0*PI/(2d0*L-1d0))
    fo(0:2*L-2,0:2*L-2) = transpose(fo(0:2*L-2,0:2*L-2)) * exp(I*(L-1)*2d0*PI/(2d0*L-1d0))
    Fmme(0:L-1,0:L-1) = fe(0:L-1, 0:L-1)
    Fmme(-(L-1):-1,0:L-1) = fe(L:2*L-2, 0:L-1)
    Fmme(0:L-1,-(L-1):-1) = fe(0:L-1, L:2*L-2)
    Fmme(-(L-1):-1,-(L-1):-1) = fe(L:2*L-2, L:2*L-2)
    Fmmo(0:L-1,0:L-1) = fo(0:L-1, 0:L-1)
    Fmmo(-(L-1):-1,0:L-1) = fo(L:2*L-2, 0:L-1)
    Fmmo(0:L-1,-(L-1):-1) = fo(0:L-1, L:2*L-2)
    Fmmo(-(L-1):-1,-(L-1):-1) = fo(L:2*L-2, L:2*L-2)

    ! Apply phase modulation to account for sampling offset.
    do m = 0, 2*L-2
       do mm = 0, 2*L-2
          Fmme(m-(L-1),mm-(L-1)) = Fmme(m-(L-1),mm-(L-1)) * exp(-I*(m+mm)*PI/(2d0*L-1d0))
          Fmmo(m-(L-1),mm-(L-1)) = Fmmo(m-(L-1),mm-(L-1)) * exp(-I*(m+mm)*PI/(2d0*L-1d0))
       end do
    end do

    Fmme(-(L-1):L-1, -(L-1):L-1) = Fmme(-(L-1):L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2
    Fmmo(-(L-1):L-1, -(L-1):L-1) = Fmmo(-(L-1):L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2

    ! Compute Gmm for even and odd extensions by direct calculation 
    ! of convolution.
    Gmme(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    Gmmo(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do m = -(L-1), L-1
       do mm = -(L-1), L-1
          do k = -(L-1), L-1 
             Gmme(m,mm) = Gmme(m,mm) + &
                  Fmme(m,k) * ssht_sampling_weight_mw(k - mm) * 2d0 * PI
             Gmmo(m,mm) = Gmmo(m,mm) + &
                  Fmmo(m,k) * ssht_sampling_weight_mw(k - mm) * 2d0 * PI
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             if (mod(m+spin,2) == 0) then
                ! m+spin even
                Gmm_term = Gmme(m,mm)
             else
                ! m+spin odd
                Gmm_term = Gmmo(m,mm)
             end if

             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Gmm_term

          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_forward_sov


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_forward_sov_conv
  !
  !! Compute forward transform for MWEO method using separation of
  !! variables, fast Fourier transforms and performing convolution with
  !! weights as product in transformed space.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_forward_sov_conv(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fe(0:2*L-2 ,0:2*L-2)
    complex(dpc) :: fo(0:2*L-2 ,0:2*L-2)
    complex(dpc) :: Fmme(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmmo(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmme(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmmo(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmm_term
    integer*8 :: fftw_plan_fwd, fftw_plan_bwd
    
    integer :: r
    complex(dpc) :: Fmme_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: Fmmo_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: w(-2*(L-1):2*(L-1))
    complex(dpc) :: wr(-2*(L-1):2*(L-1))
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_forward_sov_conv...'
       end if
    end if

    ! Extend f to the torus with even and odd extensions 
    ! about theta=PI.
    fe(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fe(L:2*L-2, 0:2*L-2) = f(L-2:0:-1, 0:2*L-2)
    fo(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fo(L:2*L-2, 0:2*L-2) = -f(L-2:0:-1, 0:2*L-2)

! If don't apply spatial shift below then apply phase modulation here.
!!$    do p = 0, 2*L-2
!!$       phi = ssht_sampling_mweo_p2phi(p, L)
!!$       do t = 0, 2*L-2
!!$          theta = ssht_sampling_mweo_t2theta(t, L)   
!!$          fe(t,p) = fe(t,p) * exp(I*(phi+theta)*(L-1))
!!$          fo(t,p) = fo(t,p) * exp(I*(phi+theta)*(L-1))
!!$       end do
!!$    end do

    ! Compute Fmm for even and odd extensions by 2D FFT.
    ! ** NOTE THAT 2D FFTW SWITCHES DIMENSIONS! HENCE TRANSPOSE BELOW. **
    call dfftw_plan_dft_2d(fftw_plan_fwd, 2*L-1, 2*L-1, fe(0:2*L-2,0:2*L-2), &
         fe(0:2*L-2,0:2*L-2), FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(fftw_plan_fwd, fe(0:2*L-2,0:2*L-2), fe(0:2*L-2,0:2*L-2))
    call dfftw_execute_dft(fftw_plan_fwd, fo(0:2*L-2,0:2*L-2), fo(0:2*L-2,0:2*L-2))
    call dfftw_destroy_plan(fftw_plan_fwd)

! If apply phase shift above just copy Fmm (and transpose 
! to account for 2D FFT).
!!$    Fmme(-(L-1):L-1, -(L-1):L-1) = transpose(fe(0:2*L-2,0:2*L-2))
!!$    Fmmo(-(L-1):L-1, -(L-1):L-1) = transpose(fo(0:2*L-2,0:2*L-2))

    ! Apply spatial shift in frequency.
    fe(0:2*L-2,0:2*L-2) = transpose(fe(0:2*L-2,0:2*L-2)) * exp(I*(L-1)*2d0*PI/(2d0*L-1d0))
    fo(0:2*L-2,0:2*L-2) = transpose(fo(0:2*L-2,0:2*L-2)) * exp(I*(L-1)*2d0*PI/(2d0*L-1d0))
    Fmme(0:L-1,0:L-1) = fe(0:L-1, 0:L-1)
    Fmme(-(L-1):-1,0:L-1) = fe(L:2*L-2, 0:L-1)
    Fmme(0:L-1,-(L-1):-1) = fe(0:L-1, L:2*L-2)
    Fmme(-(L-1):-1,-(L-1):-1) = fe(L:2*L-2, L:2*L-2)
    Fmmo(0:L-1,0:L-1) = fo(0:L-1, 0:L-1)
    Fmmo(-(L-1):-1,0:L-1) = fo(L:2*L-2, 0:L-1)
    Fmmo(0:L-1,-(L-1):-1) = fo(0:L-1, L:2*L-2)
    Fmmo(-(L-1):-1,-(L-1):-1) = fo(L:2*L-2, L:2*L-2)

    ! Apply phase modulation to account for sampling offset.
    do m = 0, 2*L-2
       do mm = 0, 2*L-2
          Fmme(m-(L-1),mm-(L-1)) = Fmme(m-(L-1),mm-(L-1)) * exp(-I*(m+mm)*PI/(2d0*L-1d0))
          Fmmo(m-(L-1),mm-(L-1)) = Fmmo(m-(L-1),mm-(L-1)) * exp(-I*(m+mm)*PI/(2d0*L-1d0))
       end do
    end do

    Fmme(-(L-1):L-1, -(L-1):L-1) = Fmme(-(L-1):L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2
    Fmmo(-(L-1):L-1, -(L-1):L-1) = Fmmo(-(L-1):L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2

    ! Compute weights.
    do mm = -2*(L-1), 2*(L-1)
       w(mm) = ssht_sampling_weight_mw(mm)
    end do

    ! Compute IFFT of w to give wr.
    wr(1:2*(L-1)) = w(-2*(L-1):-1)
    wr(-2*(L-1):0) = w(0:2*(L-1))
    w(-2*(L-1):2*(L-1)) = wr(-2*(L-1):2*(L-1))
    call dfftw_plan_dft_1d(fftw_plan_bwd, 4*L-3, wr(-2*(L-1):2*(L-1)), &
         wr(-2*(L-1):2*(L-1)), FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_execute_dft(fftw_plan_bwd, w(-2*(L-1):2*(L-1)), w(-2*(L-1):2*(L-1)))
    wr(0:2*(L-1)) = w(-2*(L-1):0)
    wr(-2*(L-1):-1) = w(1:2*(L-1))

    ! Plan forward FFT.
    call dfftw_plan_dft_1d(fftw_plan_fwd, 4*L-3, w(-2*(L-1):2*(L-1)), &
         w(-2*(L-1):2*(L-1)), FFTW_FORWARD, FFTW_MEASURE)

    ! Compute Gmm for even and odd extensions by convolution implemented 
    ! as product in real space.
    do m = -(L-1), L-1

       ! Zero-pad Fmme.
       Fmme_pad(-2*(L-1):-L) = cmplx(0d0, 0d0)
       Fmme_pad(-(L-1):L-1) = Fmme(m,-(L-1):L-1)
       Fmme_pad(L:2*(L-1)) = cmplx(0d0, 0d0)

       ! Compute IFFT of Fmme (Fmmo used for temporary storage).
       Fmmo_pad(1:2*(L-1)) = Fmme_pad(-2*(L-1):-1)
       Fmmo_pad(-2*(L-1):0) = Fmme_pad(0:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_bwd, Fmme_pad(-2*(L-1):2*(L-1)), &
            Fmme_pad(-2*(L-1):2*(L-1)))
       Fmmo_pad(0:2*(L-1)) = Fmme_pad(-2*(L-1):0)
       Fmmo_pad(-2*(L-1):-1) = Fmme_pad(1:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))

       ! Compute product of Fmme and weight in real space.
       do r = -2*(L-1), 2*(L-1)
          Fmme_pad(r) = Fmme_pad(r) * wr(-r)
       end do

       ! Compute Gmme by FFT.
       Fmmo_pad(1:2*(L-1)) = Fmme_pad(-2*(L-1):-1)
       Fmmo_pad(-2*(L-1):0) = Fmme_pad(0:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_fwd, Fmme_pad(-2*(L-1):2*(L-1)), &
            Fmme_pad(-2*(L-1):2*(L-1)))
       Fmmo_pad(0:2*(L-1)) = Fmme_pad(-2*(L-1):0)
       Fmmo_pad(-2*(L-1):-1) = Fmme_pad(1:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))

       ! Extract section of Gmme of interest.
       Gmme(m,-(L-1):(L-1)) = Fmme_pad(-(L-1):(L-1)) * 2d0 * PI / (4d0*L-3d0)

       ! Repeat for odd signal...

       ! Zero-pad Fmmo.
       Fmmo_pad(-2*(L-1):-L) = cmplx(0d0, 0d0)
       Fmmo_pad(-(L-1):L-1) = Fmmo(m,-(L-1):L-1)
       Fmmo_pad(L:2*(L-1)) = cmplx(0d0, 0d0)

       ! Compute IFFT of Fmmo (Fmme used for temporary storage).
       Fmme_pad(1:2*(L-1)) = Fmmo_pad(-2*(L-1):-1)
       Fmme_pad(-2*(L-1):0) = Fmmo_pad(0:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_bwd, Fmmo_pad(-2*(L-1):2*(L-1)), &
            Fmmo_pad(-2*(L-1):2*(L-1)))
       Fmme_pad(0:2*(L-1)) = Fmmo_pad(-2*(L-1):0)
       Fmme_pad(-2*(L-1):-1) = Fmmo_pad(1:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))

       ! Compute product of Fmmo and weight in real space.
       do r = -2*(L-1), 2*(L-1)
          Fmmo_pad(r) = Fmmo_pad(r) * wr(-r)
       end do

       ! Compute Gmmo by FFT.
       Fmme_pad(1:2*(L-1)) = Fmmo_pad(-2*(L-1):-1)
       Fmme_pad(-2*(L-1):0) = Fmmo_pad(0:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_fwd, Fmmo_pad(-2*(L-1):2*(L-1)), &
            Fmmo_pad(-2*(L-1):2*(L-1)))
       Fmme_pad(0:2*(L-1)) = Fmmo_pad(-2*(L-1):0)
       Fmme_pad(-2*(L-1):-1) = Fmmo_pad(1:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))

       ! Extract section of Gmmo of interest.
       Gmmo(m,-(L-1):(L-1)) = Fmmo_pad(-(L-1):(L-1)) * 2d0 * PI / (4d0*L-3d0)

    end do
    call dfftw_destroy_plan(fftw_plan_bwd)
    call dfftw_destroy_plan(fftw_plan_fwd)    

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el
             if (mod(m+spin,2) == 0) then
                ! m+spin even
                Gmm_term = Gmme(m,mm)
             else
                ! m+spin odd
                Gmm_term = Gmmo(m,mm)
             end if

             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Gmm_term

          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_forward_sov_conv


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_forward_sov_conv_sym
  !
  !! Compute forward transform for MWEO method using separation of
  !! variables, fast Fourier transforms, performing convolution with
  !! weights as product in transformed space and exploiting all symmetries
  !! (for complex spin signal).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_forward_sov_conv_sym(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: fe(0:2*L-2 ,0:2*L-2)
    complex(dpc) :: fo(0:2*L-2 ,0:2*L-2)
    complex(dpc) :: Fmme(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmmo(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmme(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmmo(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmm_term
    integer*8 :: fftw_plan_fwd, fftw_plan_bwd
    
    integer :: r
    complex(dpc) :: Fmme_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: Fmmo_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: w(-2*(L-1):2*(L-1))
    complex(dpc) :: wr(-2*(L-1):2*(L-1))
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_forward_sov_conv_sym...'
       end if
    end if

    ! Extend f to the torus with even and odd extensions 
    ! about theta=PI.
    fe(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fe(L:2*L-2, 0:2*L-2) = f(L-2:0:-1, 0:2*L-2)
    fo(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fo(L:2*L-2, 0:2*L-2) = -f(L-2:0:-1, 0:2*L-2)

! If don't apply spatial shift below then apply phase modulation here.
!!$    do p = 0, 2*L-2
!!$       phi = ssht_sampling_mweo_p2phi(p, L)
!!$       do t = 0, 2*L-2
!!$          theta = ssht_sampling_mweo_t2theta(t, L)   
!!$          fe(t,p) = fe(t,p) * exp(I*(phi+theta)*(L-1))
!!$          fo(t,p) = fo(t,p) * exp(I*(phi+theta)*(L-1))
!!$       end do
!!$    end do

    ! Compute Fmm for even and odd extensions by 2D FFT.
    ! ** NOTE THAT 2D FFTW SWITCHES DIMENSIONS! HENCE TRANSPOSE BELOW. **
    call dfftw_plan_dft_2d(fftw_plan_fwd, 2*L-1, 2*L-1, fe(0:2*L-2,0:2*L-2), &
         fe(0:2*L-2,0:2*L-2), FFTW_FORWARD, FFTW_ESTIMATE)
    call dfftw_execute_dft(fftw_plan_fwd, fe(0:2*L-2,0:2*L-2), fe(0:2*L-2,0:2*L-2))
    call dfftw_execute_dft(fftw_plan_fwd, fo(0:2*L-2,0:2*L-2), fo(0:2*L-2,0:2*L-2))
    call dfftw_destroy_plan(fftw_plan_fwd)

! If apply phase shift above just copy Fmm (and transpose 
! to account for 2D FFT).
!!$    Fmme(-(L-1):L-1, -(L-1):L-1) = transpose(fe(0:2*L-2,0:2*L-2))
!!$    Fmmo(-(L-1):L-1, -(L-1):L-1) = transpose(fo(0:2*L-2,0:2*L-2))

    ! Apply spatial shift in frequency.
    fe(0:2*L-2,0:2*L-2) = transpose(fe(0:2*L-2,0:2*L-2)) * exp(I*(L-1)*2d0*PI/(2d0*L-1d0))
    fo(0:2*L-2,0:2*L-2) = transpose(fo(0:2*L-2,0:2*L-2)) * exp(I*(L-1)*2d0*PI/(2d0*L-1d0))
    Fmme(0:L-1,0:L-1) = fe(0:L-1, 0:L-1)
    Fmme(-(L-1):-1,0:L-1) = fe(L:2*L-2, 0:L-1)
    Fmme(0:L-1,-(L-1):-1) = fe(0:L-1, L:2*L-2)
    Fmme(-(L-1):-1,-(L-1):-1) = fe(L:2*L-2, L:2*L-2)
    Fmmo(0:L-1,0:L-1) = fo(0:L-1, 0:L-1)
    Fmmo(-(L-1):-1,0:L-1) = fo(L:2*L-2, 0:L-1)
    Fmmo(0:L-1,-(L-1):-1) = fo(0:L-1, L:2*L-2)
    Fmmo(-(L-1):-1,-(L-1):-1) = fo(L:2*L-2, L:2*L-2)

    ! Apply phase modulation to account for sampling offset.
    do m = 0, 2*L-2
       do mm = 0, 2*L-2
          Fmme(m-(L-1),mm-(L-1)) = Fmme(m-(L-1),mm-(L-1)) * exp(-I*(m+mm)*PI/(2d0*L-1d0))
          Fmmo(m-(L-1),mm-(L-1)) = Fmmo(m-(L-1),mm-(L-1)) * exp(-I*(m+mm)*PI/(2d0*L-1d0))
       end do
    end do

    Fmme(-(L-1):L-1, -(L-1):L-1) = Fmme(-(L-1):L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2
    Fmmo(-(L-1):L-1, -(L-1):L-1) = Fmmo(-(L-1):L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2

    ! Compute weights.
    do mm = -2*(L-1), 2*(L-1)
       w(mm) = ssht_sampling_weight_mw(mm)
    end do

    ! Compute IFFT of w to give wr.
    wr(1:2*(L-1)) = w(-2*(L-1):-1)
    wr(-2*(L-1):0) = w(0:2*(L-1))
    w(-2*(L-1):2*(L-1)) = wr(-2*(L-1):2*(L-1))
    call dfftw_plan_dft_1d(fftw_plan_bwd, 4*L-3, wr(-2*(L-1):2*(L-1)), &
         wr(-2*(L-1):2*(L-1)), FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_execute_dft(fftw_plan_bwd, w(-2*(L-1):2*(L-1)), w(-2*(L-1):2*(L-1)))
    wr(0:2*(L-1)) = w(-2*(L-1):0)
    wr(-2*(L-1):-1) = w(1:2*(L-1))

    ! Plan forward FFT.
    call dfftw_plan_dft_1d(fftw_plan_fwd, 4*L-3, w(-2*(L-1):2*(L-1)), &
         w(-2*(L-1):2*(L-1)), FFTW_FORWARD, FFTW_MEASURE)

    ! Compute Gmm for even and odd extensions by convolution implemented 
    ! as product in real space.
    do m = -(L-1), L-1

       ! Zero-pad Fmme.
       Fmme_pad(-2*(L-1):-L) = cmplx(0d0, 0d0)
       Fmme_pad(-(L-1):L-1) = Fmme(m,-(L-1):L-1)
       Fmme_pad(L:2*(L-1)) = cmplx(0d0, 0d0)

       ! Compute IFFT of Fmme (Fmmo used for temporary storage).
       Fmmo_pad(1:2*(L-1)) = Fmme_pad(-2*(L-1):-1)
       Fmmo_pad(-2*(L-1):0) = Fmme_pad(0:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_bwd, Fmme_pad(-2*(L-1):2*(L-1)), &
            Fmme_pad(-2*(L-1):2*(L-1)))
       Fmmo_pad(0:2*(L-1)) = Fmme_pad(-2*(L-1):0)
       Fmmo_pad(-2*(L-1):-1) = Fmme_pad(1:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))

       ! Compute product of Fmme and weight in real space.
       do r = -2*(L-1), 2*(L-1)
          Fmme_pad(r) = Fmme_pad(r) * wr(-r)
       end do

       ! Compute Gmme by FFT.
       Fmmo_pad(1:2*(L-1)) = Fmme_pad(-2*(L-1):-1)
       Fmmo_pad(-2*(L-1):0) = Fmme_pad(0:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_fwd, Fmme_pad(-2*(L-1):2*(L-1)), &
            Fmme_pad(-2*(L-1):2*(L-1)))
       Fmmo_pad(0:2*(L-1)) = Fmme_pad(-2*(L-1):0)
       Fmmo_pad(-2*(L-1):-1) = Fmme_pad(1:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))

       ! Extract section of Gmme of interest.
       Gmme(m,-(L-1):(L-1)) = Fmme_pad(-(L-1):(L-1)) * 2d0 * PI / (4d0*L-3d0)

       ! Repeat for odd signal...

       ! Zero-pad Fmmo.
       Fmmo_pad(-2*(L-1):-L) = cmplx(0d0, 0d0)
       Fmmo_pad(-(L-1):L-1) = Fmmo(m,-(L-1):L-1)
       Fmmo_pad(L:2*(L-1)) = cmplx(0d0, 0d0)

       ! Compute IFFT of Fmmo (Fmme used for temporary storage).
       Fmme_pad(1:2*(L-1)) = Fmmo_pad(-2*(L-1):-1)
       Fmme_pad(-2*(L-1):0) = Fmmo_pad(0:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_bwd, Fmmo_pad(-2*(L-1):2*(L-1)), &
            Fmmo_pad(-2*(L-1):2*(L-1)))
       Fmme_pad(0:2*(L-1)) = Fmmo_pad(-2*(L-1):0)
       Fmme_pad(-2*(L-1):-1) = Fmmo_pad(1:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))

       ! Compute product of Fmmo and weight in real space.
       do r = -2*(L-1), 2*(L-1)
          Fmmo_pad(r) = Fmmo_pad(r) * wr(-r)
       end do

       ! Compute Gmmo by FFT.
       Fmme_pad(1:2*(L-1)) = Fmmo_pad(-2*(L-1):-1)
       Fmme_pad(-2*(L-1):0) = Fmmo_pad(0:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_fwd, Fmmo_pad(-2*(L-1):2*(L-1)), &
            Fmmo_pad(-2*(L-1):2*(L-1)))
       Fmme_pad(0:2*(L-1)) = Fmmo_pad(-2*(L-1):0)
       Fmme_pad(-2*(L-1):-1) = Fmmo_pad(1:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))

       ! Extract section of Gmmo of interest.
       Gmmo(m,-(L-1):(L-1)) = Fmmo_pad(-(L-1):(L-1)) * 2d0 * PI / (4d0*L-3d0)

    end do
    call dfftw_destroy_plan(fftw_plan_bwd)
    call dfftw_destroy_plan(fftw_plan_fwd)    

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)

          flm(ind) = flm(ind) + &
               (-1)**spin * elfactor &
               * exp(I*PION2*(m+spin)) &
               * dl(0,m) * dl(0,-spin) &
               * Gmme(m,0)

          do mm = 1, el
             if (mod(m+spin,2) == 0) then
                ! m+spin even
                Gmm_term = Gmme(m,mm) + (-1)**(m+spin)*Gmme(m,-mm)
             else
                ! m+spin odd
                Gmm_term = Gmmo(m,mm) + (-1)**(m+spin)*Gmmo(m,-mm)
             end if

             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Gmm_term

          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_forward_sov_conv_sym


  !----------------------------------------------------------------------------
  ! ssht_core_mweo_forward_sov_conv_sym_real
  !
  !! Compute forward transform for MWEO method of real scalar signal
  !! using separation of variables, fast Fourier transforms,
  !! performing convolution with weights as product in transformed
  !! space and exploiting all symmetries (including additional symmetries for
  !! real signals).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mweo_forward_sov_conv_sym_real(flm, f, L, verbosity)

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    real(dp), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    real(dp) :: fe(0:2*L-2 ,0:2*L-2)
    real(dp) :: fo(0:2*L-2 ,0:2*L-2)
    complex(dpc) :: Fmme(0:L-1, -(L-1):L-1)
    complex(dpc) :: Fmmo(0:L-1, -(L-1):L-1)
    complex(dpc) :: Gmme(0:L-1, -(L-1):L-1)
    complex(dpc) :: Gmmo(0:L-1, -(L-1):L-1)
    complex(dpc) :: Gmm_term
    integer*8 :: fftw_plan_fwd, fftw_plan_bwd
    
    integer :: r
    complex(dpc) :: Fmme_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: Fmmo_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: w(-2*(L-1):2*(L-1))
    complex(dpc) :: wr(-2*(L-1):2*(L-1))

    integer :: ind_nm
    integer :: spin
    complex(dpc) :: tmp(0:L-1,0:2*L-2)
    character(len=STRING_LEN) :: format_spec

    spin = 0

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' even-odd sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',TRUE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mweo_forward_sov_conv_sym_real...'
       end if
    end if

    ! Extend f to the torus with even and odd extensions 
    ! about theta=PI.
    ! Compute Fmm for even and odd extensions by 2D FFT.
    ! Apply spatial shift in frequency.
    call dfftw_plan_dft_r2c_2d(fftw_plan_fwd, 2*L-1, 2*L-1, fe(0:2*L-2,0:2*L-2), &
         tmp(0:L-1,0:2*L-2), FFTW_ESTIMATE)

    fe(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fe(L:2*L-2, 0:2*L-2) = f(L-2:0:-1, 0:2*L-2)
    call dfftw_execute_dft_r2c(fftw_plan_fwd, &
         transpose(fe(0:2*L-2,0:2*L-2)), tmp(0:L-1,0:2*L-2))
    Fmme(0:L-1,0:L-1) = tmp(0:L-1, 0:L-1)
    Fmme(0:L-1,-(L-1):-1) = tmp(0:L-1, L:2*L-2)

    fo(0:L-1, 0:2*L-2) = f(0:L-1, 0:2*L-2)
    fo(L:2*L-2, 0:2*L-2) = -f(L-2:0:-1, 0:2*L-2)
    call dfftw_execute_dft_r2c(fftw_plan_fwd, &
         transpose(fo(0:2*L-2,0:2*L-2)), tmp(0:L-1,0:2*L-2))
    Fmmo(0:L-1,0:L-1) = tmp(0:L-1, 0:L-1)
    Fmmo(0:L-1,-(L-1):-1) = tmp(0:L-1, L:2*L-2)

    call dfftw_destroy_plan(fftw_plan_fwd)

    Fmme(0:L-1, -(L-1):L-1) = Fmme(0:L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2
    Fmmo(0:L-1, -(L-1):L-1) = Fmmo(0:L-1, -(L-1):L-1) &
         / (2d0*L-1d0)**2

    ! Apply phase modulation to account for sampling offset.
    do m = 0, L-1
       do mm = -(L-1), L-1
          Fmme(m,mm) = Fmme(m,mm) * exp(-I*(m+mm)*PI/(2d0*L-1d0))
          Fmmo(m,mm) = Fmmo(m,mm) * exp(-I*(m+mm)*PI/(2d0*L-1d0))
       end do
    end do

    ! Compute weights.
    do mm = -2*(L-1), 2*(L-1)
       w(mm) = ssht_sampling_weight_mw(mm)
    end do

    ! Compute IFFT of w to give wr.
    wr(1:2*(L-1)) = w(-2*(L-1):-1)
    wr(-2*(L-1):0) = w(0:2*(L-1))
    w(-2*(L-1):2*(L-1)) = wr(-2*(L-1):2*(L-1))
    call dfftw_plan_dft_1d(fftw_plan_bwd, 4*L-3, wr(-2*(L-1):2*(L-1)), &
         wr(-2*(L-1):2*(L-1)), FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_execute_dft(fftw_plan_bwd, w(-2*(L-1):2*(L-1)), w(-2*(L-1):2*(L-1)))
    wr(0:2*(L-1)) = w(-2*(L-1):0)
    wr(-2*(L-1):-1) = w(1:2*(L-1))

    ! Plan forward FFT.
    call dfftw_plan_dft_1d(fftw_plan_fwd, 4*L-3, w(-2*(L-1):2*(L-1)), &
         w(-2*(L-1):2*(L-1)), FFTW_FORWARD, FFTW_MEASURE)

    ! Compute Gmm for even and odd extensions by convolution implemented 
    ! as product in real space.
    do m = 0, L-1

       ! Zero-pad Fmme.
       Fmme_pad(-2*(L-1):-L) = cmplx(0d0, 0d0)
       Fmme_pad(-(L-1):L-1) = Fmme(m,-(L-1):L-1)
       Fmme_pad(L:2*(L-1)) = cmplx(0d0, 0d0)

       ! Compute IFFT of Fmme (Fmmo used for temporary storage).
       Fmmo_pad(1:2*(L-1)) = Fmme_pad(-2*(L-1):-1)
       Fmmo_pad(-2*(L-1):0) = Fmme_pad(0:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_bwd, Fmme_pad(-2*(L-1):2*(L-1)), &
            Fmme_pad(-2*(L-1):2*(L-1)))
       Fmmo_pad(0:2*(L-1)) = Fmme_pad(-2*(L-1):0)
       Fmmo_pad(-2*(L-1):-1) = Fmme_pad(1:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))

       ! Compute product of Fmme and weight in real space.
       do r = -2*(L-1), 2*(L-1)
          Fmme_pad(r) = Fmme_pad(r) * wr(-r)
       end do

       ! Compute Gmme by FFT.
       Fmmo_pad(1:2*(L-1)) = Fmme_pad(-2*(L-1):-1)
       Fmmo_pad(-2*(L-1):0) = Fmme_pad(0:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_fwd, Fmme_pad(-2*(L-1):2*(L-1)), &
            Fmme_pad(-2*(L-1):2*(L-1)))
       Fmmo_pad(0:2*(L-1)) = Fmme_pad(-2*(L-1):0)
       Fmmo_pad(-2*(L-1):-1) = Fmme_pad(1:2*(L-1))
       Fmme_pad(-2*(L-1):2*(L-1)) = Fmmo_pad(-2*(L-1):2*(L-1))

       ! Extract section of Gmme of interest.
       Gmme(m,-(L-1):(L-1)) = Fmme_pad(-(L-1):(L-1)) * 2d0 * PI / (4d0*L-3d0)

       ! Repeat for odd signal...

       ! Zero-pad Fmmo.
       Fmmo_pad(-2*(L-1):-L) = cmplx(0d0, 0d0)
       Fmmo_pad(-(L-1):L-1) = Fmmo(m,-(L-1):L-1)
       Fmmo_pad(L:2*(L-1)) = cmplx(0d0, 0d0)

       ! Compute IFFT of Fmmo (Fmme used for temporary storage).
       Fmme_pad(1:2*(L-1)) = Fmmo_pad(-2*(L-1):-1)
       Fmme_pad(-2*(L-1):0) = Fmmo_pad(0:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_bwd, Fmmo_pad(-2*(L-1):2*(L-1)), &
            Fmmo_pad(-2*(L-1):2*(L-1)))
       Fmme_pad(0:2*(L-1)) = Fmmo_pad(-2*(L-1):0)
       Fmme_pad(-2*(L-1):-1) = Fmmo_pad(1:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))

       ! Compute product of Fmmo and weight in real space.
       do r = -2*(L-1), 2*(L-1)
          Fmmo_pad(r) = Fmmo_pad(r) * wr(-r)
       end do

       ! Compute Gmmo by FFT.
       Fmme_pad(1:2*(L-1)) = Fmmo_pad(-2*(L-1):-1)
       Fmme_pad(-2*(L-1):0) = Fmmo_pad(0:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_fwd, Fmmo_pad(-2*(L-1):2*(L-1)), &
            Fmmo_pad(-2*(L-1):2*(L-1)))
       Fmme_pad(0:2*(L-1)) = Fmmo_pad(-2*(L-1):0)
       Fmme_pad(-2*(L-1):-1) = Fmmo_pad(1:2*(L-1))
       Fmmo_pad(-2*(L-1):2*(L-1)) = Fmme_pad(-2*(L-1):2*(L-1))

       ! Extract section of Gmmo of interest.
       Gmmo(m,-(L-1):(L-1)) = Fmmo_pad(-(L-1):(L-1)) * 2d0 * PI / (4d0*L-3d0)

    end do
    call dfftw_destroy_plan(fftw_plan_bwd)
    call dfftw_destroy_plan(fftw_plan_fwd)    

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = 0, el
          call ssht_sampling_elm2ind(ind, el, m)

          flm(ind) = flm(ind) + &
               (-1)**spin * elfactor &
               * exp(I*PION2*(m+spin)) &
               * dl(0,m) * dl(0,-spin) &
               * Gmme(m,0)

          do mm = 1, el
             if (mod(m+spin,2) == 0) then
                ! m+spin even
                Gmm_term = Gmme(m,mm) + (-1)**(m+spin)*Gmme(m,-mm)
             else
                ! m+spin odd
                Gmm_term = Gmmo(m,mm) + (-1)**(m+spin)*Gmmo(m,-mm)
             end if

             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Gmm_term

          end do
       end do
    end do

    ! Set flm values for negative m using conjugate symmetry.
    do el = abs(spin), L-1
       do m = 1, el
          call ssht_sampling_elm2ind(ind, el, m)
          call ssht_sampling_elm2ind(ind_nm, el, -m)
          flm(ind_nm) = (-1)**m * conjg(flm(ind))
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mweo_forward_sov_conv_sym_real


  !----------------------------------------------------------------------------
  ! MW
  !----------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! ssht_core_mw_forward_sov_direct
  !
  !! Compute forward transform for MW method using separation of
  !! variables but only discrete (not fast) Fourier transforms.
  !!
  !! Notes:
  !!   - Used for code verification purposes (extremely slow!).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_forward_sov_direct(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmt(-(L-1):L-1, 0:2*L-2)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmm(-(L-1):L-1, -(L-1):L-1) 
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mw_forward_sov_direct...'
       end if
    end if

    ! Compute Fourier transform over phi, i.e. compute Fmt.
    Fmt(-(L-1):L-1,0:2*L-2) = cmplx(0d0, 0d0)    
    do p = 0, 2*L-2
       phi = ssht_sampling_mw_p2phi(p, L)
       do t = 0, L-1
          do m = -(L-1), L-1
             Fmt(m,t) = Fmt(m,t) + &
                  f(t,p) * exp(-I*m*phi)
          end do
       end do
    end do
    Fmt(-(L-1):L-1,0:2*L-2) = Fmt(-(L-1):L-1,0:2*L-2) / (2d0*L-1d0)

    ! Extend Fmt periodically.
    do m = -(L-1), L-1
       Fmt(m, L:2*L-2) = (-1)**(m+spin) * Fmt(m, L-2:0:-1)
    end do

    ! Compute Fmm.
    Fmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do t = 0, 2*L-2
       theta = ssht_sampling_mw_t2theta(t, L)    
       do m = -(L-1), L-1
          do mm = -(L-1), L-1
             Fmm(m,mm) = Fmm(m,mm) + &
                  Fmt(m,t) * exp(-I*mm*theta)
          end do
       end do
    end do
    Fmm(-(L-1):L-1, -(L-1):L-1) = Fmm(-(L-1):L-1, -(L-1):L-1) / (2d0*L-1d0)

    ! Compute Gmm by direct convolution.
    Gmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do m = -(L-1), L-1
       do mm = -(L-1), L-1
          do k = -(L-1), L-1 
             Gmm(m,mm) = Gmm(m,mm) + &
                  Fmm(m,k) * ssht_sampling_weight_mw(k - mm) * 2d0 * PI
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el             
             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Gmm(m,mm)
          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mw_forward_sov_direct


  !----------------------------------------------------------------------------
  ! ssht_core_mw_forward_sov
  !
  !! Compute forward transform for MW method using separation of
  !! variables and fast Fourier transforms.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_forward_sov(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmt(-(L-1):L-1, 0:2*L-2)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmm(-(L-1):L-1, -(L-1):L-1) 
    integer*8 :: fftw_plan
    complex(dpc) :: tmp(0:2*L-2)
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mw_forward_sov...'
       end if
    end if

    ! Compute Fourier transform over phi, i.e. compute Fmt.
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, Fmt(-(L-1):L-1,0), &
         Fmt(-(L-1):L-1,0), FFTW_FORWARD, FFTW_MEASURE)
    do t = 0, L-1
       call dfftw_execute_dft(fftw_plan, f(t,0:2*L-2), tmp(0:2*L-2))
       Fmt(0:L-1,t) = tmp(0:L-1)
       Fmt(-(L-1):-1,t) = tmp(L:2*L-2)
    end do
!    call dfftw_destroy_plan(fftw_plan)
    Fmt(-(L-1):L-1, 0:L-1) = Fmt(-(L-1):L-1, 0:L-1) / (2d0*L-1d0)

    ! Extend Fmt periodically.
    do m = -(L-1), L-1
       Fmt(m, L:2*L-2) = (-1)**(m+spin) * Fmt(m, L-2:0:-1)
    end do

    ! Compute Fourier transform over theta, i.e. compute Fmm.
    do m = -(L-1), L-1
       call dfftw_execute_dft(fftw_plan, Fmt(m,0:2*L-2), tmp(0:2*L-2))
       Fmm(m,0:L-1) = tmp(0:L-1)
       Fmm(m,-(L-1):-1) = tmp(L:2*L-2)
    end do
    Fmm(-(L-1):L-1, -(L-1):L-1) = Fmm(-(L-1):L-1, -(L-1):L-1) / (2d0*L-1d0)
    call dfftw_destroy_plan(fftw_plan)

    ! Apply phase modulation to account for sampling offset.
    do mm = -(L-1), L-1
       Fmm(-(L-1):L-1,mm) = Fmm(-(L-1):L-1, mm) * exp(-I*mm*PI/(2d0*L - 1d0))
    end do

    ! Compute Gmm by direct convolution.
    Gmm(-(L-1):L-1, -(L-1):L-1) = cmplx(0d0, 0d0)
    do m = -(L-1), L-1
       do mm = -(L-1), L-1
          do k = -(L-1), L-1 
             Gmm(m,mm) = Gmm(m,mm) + &
                  Fmm(m,k) * ssht_sampling_weight_mw(k - mm) * 2d0 * PI
          end do
       end do
    end do

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el             
             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Gmm(m,mm)
          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mw_forward_sov


  !----------------------------------------------------------------------------
  ! ssht_core_mw_forward_sov_conv
  !
  !! Compute forward transform for MW method using separation of
  !! variables, fast Fourier transforms and performing convolution with
  !! weights as product in transformed space.
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_forward_sov_conv(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmt(-(L-1):L-1, 0:2*L-2)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmm(-(L-1):L-1, -(L-1):L-1) 
    integer*8 :: fftw_plan
    complex(dpc) :: tmp(0:2*L-2)

    integer :: r
    complex(dpc) :: Fmm_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: tmp_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: w(-2*(L-1):2*(L-1))
    complex(dpc) :: wr(-2*(L-1):2*(L-1))
    integer*8 :: fftw_plan_fwd, fftw_plan_bwd
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mw_forward_sov_conv...'
       end if
    end if

    ! Compute Fourier transform over phi, i.e. compute Fmt.
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, Fmt(-(L-1):L-1,0), &
         Fmt(-(L-1):L-1,0), FFTW_FORWARD, FFTW_MEASURE)
    do t = 0, L-1
       call dfftw_execute_dft(fftw_plan, f(t,0:2*L-2), tmp(0:2*L-2))
       Fmt(0:L-1,t) = tmp(0:L-1)
       Fmt(-(L-1):-1,t) = tmp(L:2*L-2)
    end do
    Fmt(-(L-1):L-1, 0:L-1) = Fmt(-(L-1):L-1, 0:L-1) / (2d0*L-1d0)

    ! Extend Fmt periodically.
    do m = -(L-1), L-1
       Fmt(m, L:2*L-2) = (-1)**(m+spin) * Fmt(m, L-2:0:-1)
    end do

    ! Compute Fourier transform over theta, i.e. compute Fmm.
    do m = -(L-1), L-1
       call dfftw_execute_dft(fftw_plan, Fmt(m,0:2*L-2), tmp(0:2*L-2))
       Fmm(m,0:L-1) = tmp(0:L-1)
       Fmm(m,-(L-1):-1) = tmp(L:2*L-2)
    end do
    Fmm(-(L-1):L-1, -(L-1):L-1) = Fmm(-(L-1):L-1, -(L-1):L-1) / (2d0*L-1d0)
    call dfftw_destroy_plan(fftw_plan)

    ! Apply phase modulation to account for sampling offset.
    do mm = -(L-1), L-1
       Fmm(-(L-1):L-1,mm) = Fmm(-(L-1):L-1, mm) * exp(-I*mm*PI/(2d0*L - 1d0))
    end do

    ! Compute weights.
    do mm = -2*(L-1), 2*(L-1)
       w(mm) = ssht_sampling_weight_mw(mm)
    end do

    ! Compute IFFT of w to give wr.
    wr(1:2*(L-1)) = w(-2*(L-1):-1)
    wr(-2*(L-1):0) = w(0:2*(L-1))
    w(-2*(L-1):2*(L-1)) = wr(-2*(L-1):2*(L-1))
    call dfftw_plan_dft_1d(fftw_plan_bwd, 4*L-3, wr(-2*(L-1):2*(L-1)), &
         wr(-2*(L-1):2*(L-1)), FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_execute_dft(fftw_plan_bwd, w(-2*(L-1):2*(L-1)), w(-2*(L-1):2*(L-1)))
    wr(0:2*(L-1)) = w(-2*(L-1):0)
    wr(-2*(L-1):-1) = w(1:2*(L-1))

    ! Plan forward FFT.
    call dfftw_plan_dft_1d(fftw_plan_fwd, 4*L-3, w(-2*(L-1):2*(L-1)), &
         w(-2*(L-1):2*(L-1)), FFTW_FORWARD, FFTW_MEASURE)

    ! Compute Gmm by convolution implemented as product in real space.
    do m = -(L-1), L-1

       ! Zero-pad Fmm.
       Fmm_pad(-2*(L-1):-L) = cmplx(0d0, 0d0)
       Fmm_pad(-(L-1):L-1) = Fmm(m,-(L-1):L-1)
       Fmm_pad(L:2*(L-1)) = cmplx(0d0, 0d0)
       
       ! Compute IFFT of Fmm.
       tmp_pad(1:2*(L-1)) = Fmm_pad(-2*(L-1):-1)
       tmp_pad(-2*(L-1):0) = Fmm_pad(0:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_bwd, Fmm_pad(-2*(L-1):2*(L-1)), &
            Fmm_pad(-2*(L-1):2*(L-1)))
       tmp_pad(0:2*(L-1)) = Fmm_pad(-2*(L-1):0)
       tmp_pad(-2*(L-1):-1) = Fmm_pad(1:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))

       ! Compute product of Fmm and weight in real space.
       do r = -2*(L-1), 2*(L-1)
          Fmm_pad(r) = Fmm_pad(r) * wr(-r)
       end do

       ! Compute Gmm by FFT.
       tmp_pad(1:2*(L-1)) = Fmm_pad(-2*(L-1):-1)
       tmp_pad(-2*(L-1):0) = Fmm_pad(0:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_fwd, Fmm_pad(-2*(L-1):2*(L-1)), &
            Fmm_pad(-2*(L-1):2*(L-1)))
       tmp_pad(0:2*(L-1)) = Fmm_pad(-2*(L-1):0)
       tmp_pad(-2*(L-1):-1) = Fmm_pad(1:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))

       ! Extract section of Gmm of interest.
       Gmm(m,-(L-1):(L-1)) = Fmm_pad(-(L-1):(L-1)) * 2d0 * PI / (4d0*L-3d0)

    end do
    call dfftw_destroy_plan(fftw_plan_bwd)
    call dfftw_destroy_plan(fftw_plan_fwd)   

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)
          do mm = -el, el             
             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * Gmm(m,mm)
          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mw_forward_sov_conv


  !----------------------------------------------------------------------------
  ! ssht_core_mw_forward_sov_conv_sym
  !
  !! Compute forward transform for MW method using separation of
  !! variables, fast Fourier transforms, performing convolution with
  !! weights as product in transformed space and exploiting all symmetries
  !! (for complex spin signal).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_forward_sov_conv_sym(flm, f, L, spin, verbosity)

    integer, intent(in) :: L
    integer, intent(in) :: spin
    integer, intent(in), optional :: verbosity
    complex(dpc), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Fmt(-(L-1):L-1, 0:2*L-2)
    complex(dpc) :: Fmm(-(L-1):L-1, -(L-1):L-1)
    complex(dpc) :: Gmm(-(L-1):L-1, -(L-1):L-1) 
    integer*8 :: fftw_plan
    complex(dpc) :: tmp(0:2*L-2)

    integer :: r
    complex(dpc) :: Fmm_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: tmp_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: w(-2*(L-1):2*(L-1))
    complex(dpc) :: wr(-2*(L-1):2*(L-1))
    integer*8 :: fftw_plan_fwd, fftw_plan_bwd
    character(len=STRING_LEN) :: format_spec

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',FALSE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mw_forward_sov_conv_sym...'
       end if
    end if

    ! Compute Fourier transform over phi, i.e. compute Fmt.
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, Fmt(-(L-1):L-1,0), &
         Fmt(-(L-1):L-1,0), FFTW_FORWARD, FFTW_MEASURE)
    do t = 0, L-1
       call dfftw_execute_dft(fftw_plan, f(t,0:2*L-2), tmp(0:2*L-2))
       Fmt(0:L-1,t) = tmp(0:L-1)
       Fmt(-(L-1):-1,t) = tmp(L:2*L-2)
    end do
    Fmt(-(L-1):L-1, 0:L-1) = Fmt(-(L-1):L-1, 0:L-1) / (2d0*L-1d0)

    ! Extend Fmt periodically.
    do m = -(L-1), L-1
       Fmt(m, L:2*L-2) = (-1)**(m+spin) * Fmt(m, L-2:0:-1)
    end do

    ! Compute Fourier transform over theta, i.e. compute Fmm.
    do m = -(L-1), L-1
       call dfftw_execute_dft(fftw_plan, Fmt(m,0:2*L-2), tmp(0:2*L-2))
       Fmm(m,0:L-1) = tmp(0:L-1)
       Fmm(m,-(L-1):-1) = tmp(L:2*L-2)
    end do
    Fmm(-(L-1):L-1, -(L-1):L-1) = Fmm(-(L-1):L-1, -(L-1):L-1) / (2d0*L-1d0)
    call dfftw_destroy_plan(fftw_plan)

    ! Apply phase modulation to account for sampling offset.
    do mm = -(L-1), L-1
       Fmm(-(L-1):L-1,mm) = Fmm(-(L-1):L-1, mm) * exp(-I*mm*PI/(2d0*L - 1d0))
    end do

    ! Compute weights.
    do mm = -2*(L-1), 2*(L-1)
       w(mm) = ssht_sampling_weight_mw(mm)
    end do

    ! Compute IFFT of w to give wr.
    wr(1:2*(L-1)) = w(-2*(L-1):-1)
    wr(-2*(L-1):0) = w(0:2*(L-1))
    w(-2*(L-1):2*(L-1)) = wr(-2*(L-1):2*(L-1))
    call dfftw_plan_dft_1d(fftw_plan_bwd, 4*L-3, wr(-2*(L-1):2*(L-1)), &
         wr(-2*(L-1):2*(L-1)), FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_execute_dft(fftw_plan_bwd, w(-2*(L-1):2*(L-1)), w(-2*(L-1):2*(L-1)))
    wr(0:2*(L-1)) = w(-2*(L-1):0)
    wr(-2*(L-1):-1) = w(1:2*(L-1))

    ! Plan forward FFT.
    call dfftw_plan_dft_1d(fftw_plan_fwd, 4*L-3, w(-2*(L-1):2*(L-1)), &
         w(-2*(L-1):2*(L-1)), FFTW_FORWARD, FFTW_MEASURE)

    ! Compute Gmm by convolution implemented as product in real space.
    do m = -(L-1), L-1

       ! Zero-pad Fmm.
       Fmm_pad(-2*(L-1):-L) = cmplx(0d0, 0d0)
       Fmm_pad(-(L-1):L-1) = Fmm(m,-(L-1):L-1)
       Fmm_pad(L:2*(L-1)) = cmplx(0d0, 0d0)
       
       ! Compute IFFT of Fmm.
       tmp_pad(1:2*(L-1)) = Fmm_pad(-2*(L-1):-1)
       tmp_pad(-2*(L-1):0) = Fmm_pad(0:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_bwd, Fmm_pad(-2*(L-1):2*(L-1)), &
            Fmm_pad(-2*(L-1):2*(L-1)))
       tmp_pad(0:2*(L-1)) = Fmm_pad(-2*(L-1):0)
       tmp_pad(-2*(L-1):-1) = Fmm_pad(1:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))

       ! Compute product of Fmm and weight in real space.
       do r = -2*(L-1), 2*(L-1)
          Fmm_pad(r) = Fmm_pad(r) * wr(-r)
       end do

       ! Compute Gmm by FFT.
       tmp_pad(1:2*(L-1)) = Fmm_pad(-2*(L-1):-1)
       tmp_pad(-2*(L-1):0) = Fmm_pad(0:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_fwd, Fmm_pad(-2*(L-1):2*(L-1)), &
            Fmm_pad(-2*(L-1):2*(L-1)))
       tmp_pad(0:2*(L-1)) = Fmm_pad(-2*(L-1):0)
       tmp_pad(-2*(L-1):-1) = Fmm_pad(1:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))

       ! Extract section of Gmm of interest.
       Gmm(m,-(L-1):(L-1)) = Fmm_pad(-(L-1):(L-1)) * 2d0 * PI / (4d0*L-3d0)

    end do
    call dfftw_destroy_plan(fftw_plan_bwd)
    call dfftw_destroy_plan(fftw_plan_fwd)   

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       call ssht_dl_beta_operator(dl(-el:el,-el:el), PION2, el)
       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = -el, el
          call ssht_sampling_elm2ind(ind, el, m)

          flm(ind) = flm(ind) + &
               (-1)**spin * elfactor &
               * exp(I*PION2*(m+spin)) &
               * dl(0,m) * dl(0,-spin) &
               * Gmm(m,0)

          do mm = 1, el             
             flm(ind) = flm(ind) + &
                  (-1)**spin * elfactor &
                  * exp(I*PION2*(m+spin)) &
                  * dl(mm,m) * dl(mm,-spin) &
                  * (Gmm(m,mm) + (-1)**(m+spin)*Gmm(m,-mm))
          end do
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mw_forward_sov_conv_sym


  !----------------------------------------------------------------------------
  ! ssht_core_mw_forward_sov_conv_sym_real
  !
  !! Compute forward transform for MW method of real scalar signal
  !! using separation of variables, fast Fourier transforms,
  !! performing convolution with weights as product in transformed
  !! space and exploiting all symmetries (including additional symmetries for
  !! real signals).
  !!
  !! Variables:
  !!   - flm(0:L**2-1): Harmonic coefficients [output].
  !!   - f(0:L-1, 0:2*L-2): Function on sphere [input].
  !!   - L: Harmonic band-limit [input].
  !!   - spin: Spin number [input].
  !!   - verbosity: Verbosity flag [optional input].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  subroutine ssht_core_mw_forward_sov_conv_sym_real(flm, f, L, verbosity)

    integer, intent(in) :: L
    integer, intent(in), optional :: verbosity
    real(dp), intent(in) :: f(0:L-1 ,0:2*L-2)
    complex(dpc), intent(out) :: flm(0:L**2-1)

    integer :: p, m, t, mm, el, ind, k
    real(dp) :: theta, phi
    real(dp) :: elfactor

    real(dp) :: dl(0:L-1, 0:L-1)
    complex(dpc) :: Fmt(0:L-1, 0:2*L-2)
    complex(dpc) :: Fmm(0:L-1, -(L-1):L-1)
    complex(dpc) :: Gmm(0:L-1, -(L-1):L-1) 
    integer*8 :: fftw_plan
    complex(dpc) :: tmp(0:2*L-2)

    integer :: r
    complex(dpc) :: Fmm_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: tmp_pad(-2*(L-1):2*(L-1))
    complex(dpc) :: w(-2*(L-1):2*(L-1))
    complex(dpc) :: wr(-2*(L-1):2*(L-1))
    integer*8 :: fftw_plan_fwd, fftw_plan_bwd

    real(dp) :: tmpr(0:2*L-2)
    integer :: ind_nm
    integer :: spin
    character(len=STRING_LEN) :: format_spec

    real(dp) :: signs(0:L), ssign
    real(dp) :: dl_mm_spin, dl_0_spin
    real(dp) :: sqrt_tbl(0:2*(L-1)+1)
    integer :: eltmp

    ! Perform precomputations.
    do el = 0, 2*(L-1) + 1
       sqrt_tbl(el) = dsqrt(real(el,kind=dp))
    end do
    do m = 0, L-1, 2
       signs(m)   =  1.0_dp
       signs(m+1) = -1.0_dp
    enddo

    ! Set spin to zero.
    spin = 0
    ssign = signs(spin)

    ! Print messages depending on verbosity level.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a,a)') SSHT_PROMPT, &
               'Computing forward transform using McEwen and Wiaux', &
               ' sampling with'
          write(format_spec,'(a,i20,a,i20,a)') '(a,a,i', digit(L),',a,i', &
               digit(spin),',a)'
          write(*,trim(format_spec)) SSHT_PROMPT, &
               'parameters (L,spin,reality) = (', &
               L, ',', spin, ',TRUE)...'
       end if
       if (verbosity > 1) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Using routine ssht_core_mw_forward_sov_conv_sym_real...'
       end if
    end if

    ! Compute Fourier transform over phi, i.e. compute Fmt.
    call dfftw_plan_dft_r2c_1d(fftw_plan, 2*L-1, tmpr(0:2*L-2), &
         Fmt(0:L-1,0), FFTW_MEASURE)
    do t = 0, L-1             
       call dfftw_execute_dft_r2c(fftw_plan, f(t,0:2*L-2), tmp(0:L-1))
       Fmt(0:L-1,t) = tmp(0:L-1)
    end do
    call dfftw_destroy_plan(fftw_plan)
    Fmt(0:L-1, 0:L-1) = Fmt(0:L-1, 0:L-1) / (2d0*L-1d0)

    ! Extend Fmt periodically.
!!$    do m = 0, L-1
!!$       Fmt(m, L:2*L-2) = (-1)**(m+spin) * Fmt(m, L-2:0:-1)
!!$    end do
    do t = L, 2*L-2
       Fmt(0:L-1, t) = ssign * signs(0:L-1) * Fmt(0:L-1, 2*L-2-t)
    end do

    ! Compute Fourier transform over theta, i.e. compute Fmm.
    call dfftw_plan_dft_1d(fftw_plan, 2*L-1, tmp(0:2*L-2), &
         tmp(0:2*L-2), FFTW_FORWARD, FFTW_MEASURE)
    do m = 0, L-1
       call dfftw_execute_dft(fftw_plan, Fmt(m,0:2*L-2), tmp(0:2*L-2))
       Fmm(m,0:L-1) = tmp(0:L-1)
       Fmm(m,-(L-1):-1) = tmp(L:2*L-2)
    end do
    Fmm(0:L-1, -(L-1):L-1) = Fmm(0:L-1, -(L-1):L-1) / (2d0*L-1d0)
    call dfftw_destroy_plan(fftw_plan)

    ! Apply phase modulation to account for sampling offset.
    do mm = -(L-1), L-1
       Fmm(0:L-1,mm) = Fmm(0:L-1, mm) * exp(-I*mm*PI/(2d0*L - 1d0))
    end do

    ! Compute weights.
    do mm = -2*(L-1), 2*(L-1)
       w(mm) = ssht_sampling_weight_mw(mm)
    end do

    ! Compute IFFT of w to give wr.
    wr(1:2*(L-1)) = w(-2*(L-1):-1)
    wr(-2*(L-1):0) = w(0:2*(L-1))
    w(-2*(L-1):2*(L-1)) = wr(-2*(L-1):2*(L-1))
    call dfftw_plan_dft_1d(fftw_plan_bwd, 4*L-3, wr(-2*(L-1):2*(L-1)), &
         wr(-2*(L-1):2*(L-1)), FFTW_BACKWARD, FFTW_MEASURE)
    call dfftw_execute_dft(fftw_plan_bwd, w(-2*(L-1):2*(L-1)), w(-2*(L-1):2*(L-1)))
    wr(0:2*(L-1)) = w(-2*(L-1):0)
    wr(-2*(L-1):-1) = w(1:2*(L-1))

    ! Plan forward FFT.
    call dfftw_plan_dft_1d(fftw_plan_fwd, 4*L-3, w(-2*(L-1):2*(L-1)), &
         w(-2*(L-1):2*(L-1)), FFTW_FORWARD, FFTW_MEASURE)

    ! Compute Gmm by convolution implemented as product in real space.
    do m = 0, L-1

       ! Zero-pad Fmm.
       Fmm_pad(-2*(L-1):-L) = cmplx(0d0, 0d0)
       Fmm_pad(-(L-1):L-1) = Fmm(m,-(L-1):L-1)
       Fmm_pad(L:2*(L-1)) = cmplx(0d0, 0d0)
       
       ! Compute IFFT of Fmm.
       tmp_pad(1:2*(L-1)) = Fmm_pad(-2*(L-1):-1)
       tmp_pad(-2*(L-1):0) = Fmm_pad(0:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_bwd, Fmm_pad(-2*(L-1):2*(L-1)), &
            Fmm_pad(-2*(L-1):2*(L-1)))
       tmp_pad(0:2*(L-1)) = Fmm_pad(-2*(L-1):0)
       tmp_pad(-2*(L-1):-1) = Fmm_pad(1:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))

       ! Compute product of Fmm and weight in real space.
!!$       do r = -2*(L-1), 2*(L-1)
!!$          Fmm_pad(r) = Fmm_pad(r) * wr(-r)
!!$       end do
       Fmm_pad(-2*(L-1):2*(L-1)) = Fmm_pad(-2*(L-1):2*(L-1)) * wr(2*(L-1):-2*(L-1):-1)

       ! Compute Gmm by FFT.
       tmp_pad(1:2*(L-1)) = Fmm_pad(-2*(L-1):-1)
       tmp_pad(-2*(L-1):0) = Fmm_pad(0:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))
       call dfftw_execute_dft(fftw_plan_fwd, Fmm_pad(-2*(L-1):2*(L-1)), &
            Fmm_pad(-2*(L-1):2*(L-1)))
       tmp_pad(0:2*(L-1)) = Fmm_pad(-2*(L-1):0)
       tmp_pad(-2*(L-1):-1) = Fmm_pad(1:2*(L-1))
       Fmm_pad(-2*(L-1):2*(L-1)) = tmp_pad(-2*(L-1):2*(L-1))

       ! Extract section of Gmm of interest.
       Gmm(m,-(L-1):(L-1)) = Fmm_pad(-(L-1):(L-1)) * 2d0 * PI / (4d0*L-3d0)

    end do
    call dfftw_destroy_plan(fftw_plan_bwd)
    call dfftw_destroy_plan(fftw_plan_fwd)   

    ! Compute flm.
    flm(0::L**2-1) = cmplx(0d0, 0d0)
    do el = abs(spin), L-1
       if (el /= 0 .and. el == abs(spin)) then
          ! Recurse Wigner plane from 0 up to first el, i.e. abs(spin).
          do eltmp = 0, abs(spin)
             call ssht_dl_halfpi_trapani_eighth_table(dl(0:eltmp,0:eltmp), eltmp, &
                  sqrt_tbl(0:2*eltmp+1))
          end do
       else
          call ssht_dl_halfpi_trapani_eighth_table(dl(0:el,0:el), el, &
               sqrt_tbl(0:2*el+1))
          call ssht_dl_halfpi_trapani_fill_eighth2quarter(dl(0:el,0:el), el)
       end if

       elfactor = sqrt((2d0*el+1d0)/(4d0*PI))
       do m = 0, el
          call ssht_sampling_elm2ind(ind, el, m)

          if (spin <= 0) then
             dl_0_spin = dl(0,-spin)
          else
             dl_0_spin = signs(el) * dl(0,spin)
          end if
          
          flm(ind) = flm(ind) + &
               ssign * elfactor &
               * exp(I*PION2*(m+spin)) &
!               * dl(0,m) * dl(0,-spin) &
               * dl(0,m) * dl_0_spin &
               * Gmm(m,0)

          do mm = 1, el 

             if (spin <= 0) then
                dl_mm_spin = dl(mm,-spin)
             else
                dl_mm_spin = signs(el) * signs(mm) * dl(mm,spin)
             end if
            
             flm(ind) = flm(ind) + &
                  ssign * elfactor &
                  * exp(I*PION2*(m+spin)) &
!!$                  * dl(mm,m) * dl(mm,-spin) &
                  * dl(mm,m) * dl_mm_spin &
                  * (Gmm(m,mm) + ssign * signs(m) * Gmm(m,-mm))
          end do
       end do
    end do

    ! Set flm values for negative m using conjugate symmetry.
    do el = abs(spin), L-1
       do m = 1, el
          call ssht_sampling_elm2ind(ind, el, m)
          call ssht_sampling_elm2ind(ind_nm, el, -m)
          flm(ind_nm) = signs(m) * conjg(flm(ind))
       end do
    end do

    ! Print finished if verbosity set.
    if (present(verbosity)) then
       if (verbosity > 0) then
          write(*,'(a,a)') SSHT_PROMPT, &
               'Forward transform computed!'
       end if
    end if

  end subroutine ssht_core_mw_forward_sov_conv_sym_real



  !--------------------------------------------------------------------------
  ! Utility routines
  !--------------------------------------------------------------------------


  !----------------------------------------------------------------------------
  ! digit
  !
  !! Compute number of digits required to display integer (must be
  !! less than 25).
  !!
  !! Variables:
  !!   - L: Integer considered [input].
  !!   - d: Number of digits [output].
  !
  !! @author J. D. McEwen
  !
  ! Revisions:
  !   October 2010 - Written by Jason McEwen
  !----------------------------------------------------------------------------

  function digit(L) result(d)

    integer, intent(in) :: L
    integer :: d

    if(L/10 < 1) then
       d = 1
    elseif(L/100 < 1) then
       d = 2
    elseif(L/1000 < 1) then
       d = 3
    elseif(L/10000 < 1) then
       d = 4
    elseif(L/100000 < 1) then
       d = 5
    else
       d = 25
    end if
 
  end function digit


end module ssht_core_mod
