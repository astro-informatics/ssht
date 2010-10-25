!------------------------------------------------------------------------------
! ssht_wav2sky
!
!! Converts wavelet coefficients read from a SSHT formatted fits/matlab file 
!! (equi-angular pixelisation) to a sky Healpix fits file.  Either the sky corresponding to 
!! a specified dilation and orientation index may be written only, or skies 
!! corresponding to all dilations and orientations may be written.
!!
!! Usage: ssht_wav2sky
!!   - [-help]: Displays usage information.
!!   - [-inp filename_in]: Name of input wavelet coefficient SSHT formatted 
!!     fits/matlab file.
!!   - [-file_type file_type (fits; m)]: String specifying type of input 
!!     SSHT file to read  (fits or matlab m file) [default=fits].
!!   - [-out filename_out]: Name of output sky coefficient map (only required 
!!     if all=.false.)
!!   - [-jj jj]: Scale index of coefficients to write (only 
!!     required if all=.false.)
!!   - [-gg gg]: Orientation index of coefficients to write
!!     (only required if all=.false.)
!!   - [-nside healpix_nside]: Optional Healpix nside resolution of written sky. 
!!     If this option is *not* specified then an appropriate nside value is 
!!     set for each level jj.
!!   - [-interp interpolation_status]: Logical to specify whether to perform
!!     linear interpolation (if true) or else simply take nearest pixel
!!     (default=false).
!!   - [-all all_status]: Logical to specify whether skies corresponding to 
!!     all dilations and orientations are to be written (default=false).
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2007
!
! Revisions:
!   November 2007 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program ssht_wav2sky

	use ssht_types_mod
	use ssht_error_mod
	use ssht_fileio_mod
	use ssht_core_mod
	use s2_sky_mod

  implicit none

  type(ssht_wav_abg), allocatable :: wavdyn(:)
	complex(dpc), allocatable :: scoeff(:,:)
	integer :: J
	integer :: B
	integer :: N
	integer :: bl_scoeff
	real(dp) :: alpha
	character(len=STRING_LEN) :: filename_in, filename_out

	type(s2_sky) :: sky
	integer :: jj, gg, bl_hi, nside, pix_scheme
	logical :: interp, all
	logical :: use_default_nside
  integer :: SSHT_FITS_FILENAME_EXT_LEN = 4

  character(len=*), parameter ::  FILE_TYPE_FITS = 'fits'
  character(len=*), parameter ::  FILE_TYPE_MAT = 'm'
  character(len=STRING_LEN) :: file_type = FILE_TYPE_FITS

	! Set default parameter values.
	filename_out = 'sky_out.fits'
  nside = 128
  interp = .false.
  all = .false.
	use_default_nside = .true.
	pix_scheme = S2_SKY_RING
	jj = 0
	gg = 0

	! Parse options from command line.
	call parse_options()

	! Read SSHT fits file.
  select case (trim(file_type))

    case (FILE_TYPE_FITS)
       call ssht_fileio_fits_wav_read(wavdyn, scoeff, J, B, N, bl_scoeff, &
            alpha, filename_in)
       SSHT_FITS_FILENAME_EXT_LEN = 5

    case (FILE_TYPE_MAT)
       call ssht_fileio_matlab_wav_read(wavdyn, scoeff, J, B, N, &
            bl_scoeff, alpha, filename_in)
       SSHT_FITS_FILENAME_EXT_LEN = 2

    case default
       call ssht_error(SSHT_ERROR_FILEIO, 'ssht_wav2sky', &
            comment_add='Invalid file type option')

  end select

	if(all) then

     ! Write all maps.
     do jj = 0,J
       do gg = 0,N-1

         if(J < 100) then
           write(filename_out, '(a,a,i2.2,a,i2.2,a)') &
             trim(filename_in(1:len(trim(filename_in)) &
                              - SSHT_FITS_FILENAME_EXT_LEN)), &
             '_sky_jj', jj, '_gg', gg, '.fits'
          else
            write(filename_out, '(a,a,i3.3,a,i2.2,a)') &
              trim(filename_in(1:len(trim(filename_in)) &
                                - SSHT_FITS_FILENAME_EXT_LEN)), &
              '_sky_jj', jj, '_gg', gg, '.fits'
         end if

				! Compute band limit for given scale.
				bl_hi = min(ceiling(B / (alpha**(jj-1)) - TOL_CEIL) , B)

				!if(use_default_nside) nside = bl_hi/2

				! Create sky from wavelet coefficients.
				sky = s2_sky_init(real(wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, gg),sp), interp, &
					nside, pix_scheme, sdw=.true.)

				! Write sky to file.
				call s2_sky_write_map_file(sky, filename_out, &
					comment='Display map computed from SSHT')
		
				! Free temporary sky extracted.
				call s2_sky_free(sky)

        end do
     end do

	else

		! Check jj and gg valid.
		if(gg<0 .or. gg>=N .or. jj<0 .or. jj>J) then
			call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_wav2sky', &
				comment_add='Scale or orientation index out of bounds')
		end if

		! Compute band limit for given scale.
		bl_hi = min(ceiling(B / (alpha**(jj-1)) - TOL_CEIL) , B)

		!if(use_default_nside) nside = bl_hi/2

		! Create sky from wavelet coefficients.
		sky = s2_sky_init(real(wavdyn(jj)%coeff(0:2*bl_hi-2, 0:2*bl_hi-1, gg),sp), interp, &
			nside, pix_scheme, sdw=.true.)

		! Write sky to file.
		call s2_sky_write_map_file(sky, filename_out, &
			comment='Display map computed from SSHT')

		! Free temporary sky extracted.
		call s2_sky_free(sky)

	end if

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
      
      integer :: n, i
      character(len=STRING_LEN) :: opt
      character(len=STRING_LEN) :: arg
      
      n = nArguments()
     
      do i = 1,n,2
        
        call getArgument(i,opt)
     
        if (i == n .and. trim(opt) /= '-help') then
          write(*,*) 'Error: Option ', trim(opt), ' has no argument'
          stop
        end if
     
        if(trim(opt) /= '-help') call getArgument(i+1,arg)

        ! Read each argument in turn
        select case (trim(opt))
  
          case ('-help')
            write(*,'(a)') 'Usage: ssht_wav2sky [-inp filename_in]'
            write(*,'(a)') '                   [-file_type file_type (fits; m)]'
            write(*,'(a)') '                   [-out filename_out]'
            write(*,'(a)') '                   [-all all_status]'  
            write(*,'(a)') '                   [-jj jj in range [0,J]]' 
            write(*,'(a)') '                   [-gg gg in range [0,N-1]]' 
            write(*,'(a)') '                   [-nside healpix_nside]' 
            write(*,'(a)') &
              '                   [-interp interpolation_status (optional)]' 

            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-file_type')
            file_type = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-jj')
            read(arg,*) jj

          case ('-gg')
            read(arg,*) gg

         case ('-nside')
            read(arg,*) nside
						!use_default_nside = .false.

         case ('-interp')
            read(arg,*) interp

          case ('-all')
            read(arg,*) all

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program ssht_wav2sky

