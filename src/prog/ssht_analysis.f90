!------------------------------------------------------------------------------
! ssht_analysis
!
!! Computes the SSHT wavelet and scaling coefficients of a Healpix sky map.
!!
!! Notes:
!!   - The harmonic band limit B is set to 2*nside, where nside is the 
!!     Healpix resolution parameter of the input sky map.
!!   - The spherical harmonic transform and inverse provided by Healpix 
!!     is not exact, hence these reduce the numerical precision of the 
!!     `exact' reconstruction of the original real space map from its wavelet 
!!     coefficients.
!! 
!! Usage: ssht_analysis
!!   - [-help]: Displays usage information.
!!   - [-inp filename_in]: Name of input Healpix sky fits map.
!!   - [-out filename_out]: Name of output SSHT formatted fits/matlab file 
!!     containing wavelet and scaling coefficients.
!!   - [-file_type file_type (fits; m)]: String specifying type of output 
!!     SSHT file to write  (fits or matlab m) [default=fits].
!!   - [-N N]: Azimuthal band limit [default=3].
!!   - [-alpha alpha]: Basis dilation factor [default=2].
!!   - [-J J]: Maximum analysis scale (optional) [default=Jmax].
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2007
!
! Revisions:
!   November 2007 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program ssht_analysis

	use ssht_types_mod
	use ssht_error_mod
	use ssht_fileio_mod
	use ssht_core_mod
	use s2_sky_mod

  implicit none

	complex(dpc), allocatable :: flm(:,:)
	complex(spc), allocatable :: flm_temp(:,:)
	real(dp), allocatable :: K_gamma(:,:)
	real(dp), allocatable :: Phi2(:)
	complex(dpc), allocatable :: Slm(:,:)
	real(dp), allocatable :: admiss(:)
  type(ssht_wav_abg), allocatable :: wavdyn(:)
	complex(dpc), allocatable :: scoeff(:,:)
	integer :: J
	integer :: J_max
	integer :: B
	integer :: N
	integer :: bl_scoeff
	real(dp) :: alpha
	integer :: nside
	logical :: admiss_pass
	type(s2_sky) :: sky
	integer :: fail = 0
	logical :: use_Jmax
	character(len=STRING_LEN) :: filename_in, filename_out
  integer :: file_extension = 1

  character(len=*), parameter ::  FILE_TYPE_FITS = 'fits'
  character(len=*), parameter ::  FILE_TYPE_MAT = 'm'
  character(len=STRING_LEN) :: file_type = FILE_TYPE_FITS
  character(len=STRING_LEN) :: error_string

	! Set default parameter values.
	N = 3
	alpha = 2d0
	use_Jmax = .true.
	filename_in = 'wmap_ilc_3yr_v2_n64.fits'
	filename_out = 'wmap_ilc_3yr_v2_n64.ssht'

	! Parse options from command line.
	call parse_options()

	! Read sky.
	sky = s2_sky_init(filename_in, file_extension)
	nside = s2_sky_get_nside(sky)
	B = 2*nside
	J_max = ceiling(log(real(B,dp))/log(alpha) - TOL_CEIL)
	if(use_Jmax) J = J_max
  if(J > J_max) then
     J = J_max
     write(error_string,'(a,i4)') 'J too large, setting to maximum J = ', J_max
     call ssht_error(SSHT_ERROR_ARG_WARNING, 'ssht_analysis', &
          comment_add=trim(error_string))
  end if

	! Allocate memory.
	allocate(flm_temp(0:B-1,0:B-1), stat=fail)
	allocate(flm(0:B-1,0:B-1), stat=fail)
	allocate(K_gamma(0:J,0:B-1), stat=fail)
	allocate(Phi2(0:B-1), stat=fail)
	allocate(Slm(0:B-1,0:N-1), stat=fail)
	allocate(admiss(0:B-1), stat=fail)
	if(fail /= 0) then
		call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_analysis')
	end if

	! Compute spherical harmonic coefficients.
	call s2_sky_compute_alm(sky, B-1, B-1)
	call s2_sky_get_alm(sky, flm_temp(0:B-1,0:B-1))
	flm(0:B-1,0:B-1) = flm_temp(0:B-1,0:B-1)

	! Compute kernels, scaling function and directionality coefficients.
	call ssht_core_init_kernels(K_gamma, Phi2, bl_scoeff, J, B, alpha)
	call ssht_core_init_directionality(Slm, B, N)
	admiss_pass = ssht_core_admiss(admiss, K_gamma, Phi2, B, J)
	if(.not. admiss_pass) then
		call ssht_error(SSHT_ERROR_ADMISS_FAIL, 'ssht_analysis')
	end if

	! Allocate memory for scaling coefficients (cannot be done earlier 
	! since don't know bl_scoeff).
	allocate(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
	if(fail /= 0) then
		call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_analysis')
	end if

	! Perform SSHT analysis.
	call ssht_core_analysis_flm2wav_dynamic(wavdyn, scoeff, flm, K_gamma, Slm, &
		J, B, N, bl_scoeff, alpha)

	! Save wavelet and scaling coefficients.
  select case (trim(file_type))

    case (FILE_TYPE_FITS)
       call ssht_fileio_fits_wav_write(wavdyn, scoeff, J, B, N, bl_scoeff, &
            alpha, filename_out)

    case (FILE_TYPE_MAT)
       call ssht_fileio_matlab_wav_write(wavdyn, scoeff, J, B, N, &
            bl_scoeff, alpha, filename_out)

    case default
       call ssht_error(SSHT_ERROR_FILEIO, 'ssht_analysis', &
            comment_add='Invalid file type option')

  end select

	! Free memory.
	deallocate(flm_temp, flm, K_gamma, Phi2, Slm, admiss, scoeff)
  call ssht_core_free_wavdyn(wavdyn)
	call s2_sky_free(sky)

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
            write(*,'(a)') 'Usage: ssht_analysis [-inp filename_in]'
            write(*,'(a)') '                    [-out filename_out]'
            write(*,'(a)') '                    [-file_type file_type (fits; m)]'
            write(*,'(a)') '                    [-N N]'  
            write(*,'(a)') '                    [-alpha alpha]' 
            write(*,'(a)') '                    [-J J]'
            stop
          
          case ('-inp')
            filename_in = trim(arg)

          case ('-out')
            filename_out = trim(arg)

          case ('-file_type')
            file_type = trim(arg)

          case ('-N')
            read(arg,*) N

          case ('-alpha')
            read(arg,*) alpha

         case ('-J')
            read(arg,*) J
						use_Jmax = .false.

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program ssht_analysis
