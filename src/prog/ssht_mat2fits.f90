!------------------------------------------------------------------------------
! ssht_mat2fits
!
!! Converts a matlab SSHT file containing wavelet and scaling coefficients 
!! to a fits SSHT file.
!! 
!! Usage: ssht_mat2fits
!!   - [-help]: Displays usage information.
!!   - [-inp filename_in]: Name of input matlab SSHT m file.
!!   - [-out filename_out]: Name of output fits SSHT file.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - May 2008
!
! Revisions:
!   May 2008 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program ssht_mat2fits

	use ssht_types_mod
	use ssht_error_mod
	use ssht_fileio_mod
  use ssht_core_mod

  implicit none

  type(ssht_wav_abg), allocatable :: wavdyn(:)
	complex(dpc), allocatable :: scoeff(:,:)
	integer :: J
	integer :: B
	integer :: N
	integer :: bl_scoeff
	real(dp) :: alpha
	character(len=STRING_LEN) :: filename_in, filename_out

	! Parse options from command line.
	call parse_options()

  ! Read wavelet and scaling coefficients from matlab file.
  call ssht_fileio_matlab_wav_read(wavdyn, scoeff, J, B, N, &
       bl_scoeff, alpha, filename_in)

  ! Write wavelet and scaling coefficients to fits file.
  call ssht_fileio_fits_wav_write(wavdyn, scoeff, J, B, N, bl_scoeff, &
       alpha, filename_out)

	! Free memory.
	deallocate(scoeff)
  call ssht_core_free_wavdyn(wavdyn)

	!----------------------------------------------------------------------------

  contains


    !---------------------------------------------------------------------
    ! parse_options
    !
    !! Parses the options passed when program called.
    !
    !!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
    !! @version 0.1 - November 2007
    !
    ! Revisions:
    !   November 2007 - Written by Jason McEwen 
    !---------------------------------------------------------------------

    subroutine parse_options()

      use extension, only: getArgument, nArguments
     
      implicit none
      
      integer :: nn, i
      character(len=STRING_LEN) :: opt
      character(len=STRING_LEN) :: arg
      
      nn = nArguments()
     
      do i = 1,nn,2
        
        call getArgument(i,opt)
     
        if (i == nn .and. trim(opt) /= '-help') then
          write(*,*) 'Error: Option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: ssht_mat2fits [-inp filename_in]'
            write(*,'(a)') '                    [-out filename_out]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program ssht_mat2fits
