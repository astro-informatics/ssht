!------------------------------------------------------------------------------
! ssht_fileio_mod  -- SSHT library file IO class
! 
!! Functionality to read and write SSHT formatted fits and matlab files 
!! containing wavelet and scaling coefficients.  Both statically and 
!! dynamically allocated wavelet coefficients may be written and read from 
!! files (both data types have the same fits file format).
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 November 2007
!
! Revisions:
!   November 2007 - Written by Jason McEwen
!------------------------------------------------------------------------------

module ssht_fileio_mod

  use ssht_types_mod
  use ssht_error_mod
  use ssht_core_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  integer :: FILEIO = 1

  !---------------------------------------
  ! Interfaces
  !---------------------------------------



  !----------------------------------------------------------------------------

contains


  !--------------------------------------------------------------------------
  ! Matlab file IO 
  !--------------------------------------------------------------------------



end module ssht_fileio_mod



