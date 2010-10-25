!------------------------------------------------------------------------------
! ssht_wavplot
!
!! Computes a Healpix sky map of wavelet for a given j for subsequent 
!! plotting.
!! 
!! Usage: ssht_wavplot
!!   - [-help]: Displays usage information.
!!   - [-nside nside]: Healpix resolution parameter [default=256].
!!   - [-jj jj]: Analysis depth to plot wavelet at [default=0].
!!   - [-alpha alpha]: Basis dilation factor [default=2].
!!   - [-B B]: Harmonic band lmit [default=128].
!!   - [-N N]: Azimuthal band limit [default=3].
!!   - [-out filename_out]: Name of output Healpix fits file
!!     containing wavelet.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2007
!
! Revisions:
!   November 2007 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program ssht_wavplot

	use ssht_types_mod
	use ssht_error_mod
	use ssht_core_mod
	use s2_sky_mod

  implicit none

	complex(dpc), allocatable :: flm(:,:)
	complex(spc), allocatable :: flm_temp(:,:)
	real(dp), allocatable :: K_gamma(:,:)
	real(dp), allocatable :: Phi2(:)
	complex(dpc), allocatable :: Slm(:,:)
	real(dp), allocatable :: admiss(:)
	integer :: J
	integer :: J_max
	integer :: B
	integer :: N
	integer :: bl_scoeff
	real(dp) :: alpha
	integer :: nside
	type(s2_sky) :: sky
	integer :: fail = 0, el, m, jj
	character(len=STRING_LEN) :: filename_out
  real(sp), parameter :: ALPHA_CENTER = 0.0e0
  real(sp), parameter :: BETA_CENTER = pi/2.0e0
  real(sp), parameter :: GAMMA_CENTER = pi
	logical :: rotate, admiss_pass

	jj=0
	nside = 256
	B = 128
	alpha = 2d0
	N = 3
	filename_out = 'wav.fits'
	rotate = .true.

	! Parse options from command line.
	call parse_options()
 
	J_max = ssht_core_comp_Jmax(B, alpha)
	J = J_max

	! Check jj valid.
	if(jj>J_max .or. jj<0) then
		call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_wavplot', &
			comment_add='Analysis depth j invalid')
	end if

	! Allocate memory.
	allocate(flm_temp(0:B-1,0:B-1), stat=fail)
	allocate(flm(0:B-1,0:B-1), stat=fail)
	allocate(K_gamma(0:J,0:B-1), stat=fail)
	allocate(Phi2(0:B-1), stat=fail)
	allocate(Slm(0:B-1,0:N-1), stat=fail)
	allocate(admiss(0:B-1), stat=fail)
	if(fail /= 0) then
		call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_wavplot')
	end if

	! Compute kernels, scaling function and directionality coefficients.
	call ssht_core_init_kernels(K_gamma, Phi2, bl_scoeff, J, B, alpha)
	call ssht_core_init_directionality(Slm, B, N)
	admiss_pass = ssht_core_admiss(admiss, K_gamma, Phi2, B, J)
	if(.not. admiss_pass) then
		call ssht_error(SSHT_ERROR_ADMISS_FAIL, 'ssht_wavplot')
	end if

	! Compute wavelet in real space.
	do el=0,B-1
		do m = 0,min(N-1,el)
			flm(el,m) = K_gamma(jj,el) * Slm(el,m)
		end do
	end do
	flm_temp(0:B-1,0:B-1) = flm(0:B-1,0:B-1)
	sky = s2_sky_init(flm_temp(0:B-1,0:B-1), B-1, B-1)
	call s2_sky_compute_map(sky, nside)

	! Rotate to equator if required.
	if(rotate) call s2_sky_rotate(sky, ALPHA_CENTER, BETA_CENTER, GAMMA_CENTEr)

	! Save wavelet fits file.
	call s2_sky_write_map_file(sky, filename_out)
	
	! Free memory.
	call s2_sky_free(sky)
	deallocate(flm, flm_temp, K_gamma, Phi2, Slm, admiss)

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
            write(*,'(a)') 'Usage: ssht_wavplot [-nside nside]'
            write(*,'(a)') '                   [-jj jj]'
            write(*,'(a)') '                   [-alpha alpha]'
            write(*,'(a)') '                   [-B B]'
            write(*,'(a)') '                   [-N N]'
            write(*,'(a)') '                   [-out filename_out]'
            stop
          
          case ('-nside')
            read(arg,*) nside

          case ('-jj')
            read(arg,*) jj

          case ('-alpha')
            read(arg,*) alpha

          case ('-B')
            read(arg,*) B

          case ('-N')
            read(arg,*) N

          case ('-out')
            filename_out = trim(arg)

          case default
            print '("unknown option ",a4," ignored")', opt            

        end select
      end do

    end subroutine parse_options


end program ssht_wavplot
