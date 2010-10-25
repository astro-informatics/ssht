!------------------------------------------------------------------------------
! ssht_test
!
!! Performs SSHT transform analysis and synthesis and check that the original 
!! signal is reconstructed exactly (to numerical precision).  Test is 
!! performed on a random signal with harmonic coefficients uniformly 
!! sampled from (-1,1).
!!
!! Usage: ssht_test B, e.g. ssht_test 64
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2007
!
! Revisions:
!   November 2007 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program ssht_test

  use ssht_types_mod
  use ssht_error_mod
  use ssht_core_mod
  use ssht_fileio_mod
  !use F90_UNIX_ENV

  implicit none

	 interface
     subroutine ssht_test_gen_flm(L, flm, seed)
       use ssht_types_mod, only: dpc
       implicit none
       integer, intent(in) :: L
       complex(dpc), intent(out) :: flm(0:L,0:L)
       integer, intent(in) :: seed
     end subroutine ssht_test_gen_flm
  end interface

	character(len=64) :: arg
	integer, parameter :: N_repeat = 1
	integer :: B
	integer :: N
	real(dp) :: alpha
	integer :: J
	integer :: J_max
	integer :: bl_scoeff
	integer :: fail = 0, seed, i_repeat
	real(dp) :: error_flm(0:N_repeat-1)
	logical :: admiss_pass
	real :: time_start, time_end
	real :: durations_analysis(0:N_repeat-1)
	real :: durations_synthesis(0:N_repeat-1)

	complex(dpc), allocatable :: flm_orig(:,:), flm_syn(:,:)
	real(dp), allocatable :: K_gamma(:,:)
	real(dp), allocatable :: Phi2(:)
	complex(dpc), allocatable :: Slm(:,:)
	real(dp), allocatable :: admiss(:)
  type(ssht_wav_abg), allocatable :: wavdyn(:)
!  real(dp), allocatable :: wav(:,:,:,:)
	complex(dpc), allocatable :: scoeff(:,:)

	! Initialise parameters.
	call getarg(1, arg)
	read(arg,*) B
	!B = 32
	N = 3
	alpha = 2d0	
	J_max = ssht_core_comp_Jmax(B, alpha)
	J = J_max
	seed = 1

	! Allocate memory.
	allocate(flm_orig(0:B-1, 0:B-1), stat=fail)
	allocate(flm_syn(0:B-1, 0:B-1), stat=fail)
	allocate(K_gamma(0:J,0:B-1), stat=fail)
	allocate(Phi2(0:B-1), stat=fail)
	allocate(Slm(0:B-1,0:N-1), stat=fail)
	allocate(admiss(0:B-1), stat=fail)
!	allocate(wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1), stat=fail)
	if(fail /= 0) then
		call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_test')
	end if

	write(*,*)
	write(*,'(a)') 'SSHT test program'
	write(*,'(a)') '================'
	write(*,*)

	do i_repeat = 0,N_repeat-1

		! Generate harmonic coefficients of random test signal.
		call ssht_test_gen_flm(B-1, flm_orig, seed)
	
		! Compute kernels, scaling function and directionality coefficients.
		write(*,'(a,i2)') 'Initialisation no.', i_repeat
		call cpu_time(time_start)
		call ssht_core_init_kernels(K_gamma, Phi2, bl_scoeff, J, B, alpha)
		call ssht_core_init_directionality(Slm, B, N)
		admiss_pass = ssht_core_admiss(admiss, K_gamma, Phi2, B, J)
		if(.not. admiss_pass) then
			call ssht_error(SSHT_ERROR_ADMISS_FAIL, 'ssht_test')
		end if
		call cpu_time(time_end)
		write(*,'(a,f43.2)') ' duration =', time_end - time_start

		! Allocate memory for scaling coefficients (cannot be done earlier 
		! since don't know bl_scoeff).
		allocate(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
		if(fail /= 0) then
			call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_test')
		end if
	
		! Compute wavelet and scaling coefficients.
		write(*,'(a,i2)') 'Analysis no.', i_repeat
		call cpu_time(time_start)
		call ssht_core_analysis_flm2wav_dynamic(wavdyn, scoeff, flm_orig, K_gamma, Slm, &
			J, B, N, bl_scoeff, alpha)
!		call ssht_core_analysis_flm2wav(wav, scoeff, flm_orig, K_gamma, Slm, &
!			J, B, N, bl_scoeff, alpha)
		call cpu_time(time_end)
		durations_analysis(i_repeat) = time_end - time_start
		write(*,'(a,f43.2)') ' duration =', durations_analysis(i_repeat)

		! Synthesis harmonic coefficients of signal from wavelet and scaling 
		! coefficients.
		write(*,'(a,i2)') 'Synthesis no.', i_repeat
		call cpu_time(time_start)
		call ssht_core_synthesis_wav2flm_dynamic(flm_syn, wavdyn, scoeff, K_gamma, Phi2, &
				Slm, J, B, N, bl_scoeff, alpha)
!		call ssht_core_synthesis_wav2flm(flm_syn, wav, scoeff, K_gamma, Phi2, &
!				Slm, J, B, N, bl_scoeff, alpha)
		call cpu_time(time_end)
		durations_synthesis(i_repeat) = time_end - time_start
		write(*,'(a,f43.2)') ' duration =', durations_synthesis(i_repeat)
	
		error_flm(i_repeat) = maxval(abs(flm_orig(0:B-1,0:B-1) - flm_syn(0:B-1,0:B-1)))
		write(*,'(a,e43.5)') ' error =   ', error_flm(i_repeat)
		write(*,*)

		deallocate(scoeff)

	end do

	write(*,'(a)') 'Summary:'
	write(*,'(a,i30)') 'N_repeat               =', N_repeat
	write(*,'(a,i30)') 'B                      =', B
	write(*,'(a,i30)') 'N                      =', N
	write(*,'(a,i30)') 'J                      =', J	
	write(*,'(a,f30.5)') 'alpha                  =', alpha
	write(*,'(a,f30.5)') 'Average analysis time  =', sum(durations_analysis(0:N_repeat-1)) / real(N_repeat)
	write(*,'(a,f30.5)') 'Average synthesis time =', sum(durations_synthesis(0:N_repeat-1)) / real(N_repeat)
	write(*,'(a,e30.5)') 'Average max error      =', sum(error_flm(0:N_repeat-1)) / real(N_repeat)
	write(*,*)

	if( abs(sum(error_flm(0:N_repeat-1)) / real(N_repeat)) < 1d-6) then
		write(*,'(a)') 'Tests passed!'
	else
		write(*,'(a)') 'Tests failed!'
	end if 
	write(*,*) 

	! Deallocate memory.
	deallocate(flm_orig, flm_syn, K_gamma, Phi2, Slm, admiss)
!  deallocate(wav)

end program ssht_test


!--------------------------------------------------------------------------
! ssht_test_gen_flm
!
!! Generate random spherical harmonic coefficients.
!!
!! Notes: 
!!   - Uniform deviate (Num rec 1992, chap 7.1), original routine
!!     said to be 'perfect'.
!!
!! Variables:
!!   - L: Maximum harmonic and limit to consider (i.e. B-1) [input].
!!   - flm(0:L,0:L): Random spherical harmonic coefficients generated 
!!     [output].
!!   - seed: Integer seed required for random number generator [input].
!
!! @author J. D. McEwen
!! @version 0.1 October 2007
! 
! Revisions:
!   October 2007 - Jason McEwen
!--------------------------------------------------------------------------

subroutine ssht_test_gen_flm(L, flm, seed)

  use ssht_types_mod, only: dpc

  implicit none

  interface 
     function ran2_dp(idum)
       use ssht_types_mod, only: dp
       real(dp) :: ran2_dp
       integer :: idum
     end function ran2_dp
  end interface

  integer, intent(in) :: L
  complex(dpc), intent(out) :: flm(0:L,0:L)
  integer, intent(in) :: seed

  integer :: el, m

  flm(0:L,0:L) = 0d0

  do el = 0,L

     flm(el,0) = cmplx(2d0*ran2_dp(seed)-1d0, 0d0)

     do m = 1,el
        flm(el,m) = cmplx(2d0*ran2_dp(seed)-1d0, 2d0*ran2_dp(seed)-1d0)
     end do

  end do

end subroutine ssht_test_gen_flm


!--------------------------------------------------------------------------
! ran2_dp
!
!! Generate uniform deviate in range [0,1) given seed.
!! (Using double precision.)
!!
!! Notes: 
!!   - Uniform deviate (Num rec 1992, chap 7.1), original routine
!!     said to be 'perfect'.
!!
!! Variables:
!!   - idum: Seed.
!
!! @author J. D. McEwen
!! @version 0.1 May 2007
! 
! Revisions:
!   May 2007 - Integrated by Jason McEwen
!--------------------------------------------------------------------------

function ran2_dp(idum)

  use ssht_types_mod, only: dp

	implicit none

  INTEGER :: idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
  REAL(dp) :: ran2_dp,AM,EPS,RNMX
  PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
       & IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
       & NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
  INTEGER :: idum2,j,k,iv(NTAB),iy
  SAVE iv,iy,idum2
  DATA idum2/123456789/, iv/NTAB*0/, iy/0/
  if (idum.le.0) then
     idum=max(-idum,1)
     idum2=idum
     do j=NTAB+8,1,-1
        k=idum/IQ1
        idum=IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        if (j.le.NTAB) iv(j)=idum
     end do
     iy=iv(1)
  end if
  k=idum/IQ1
  idum=IA1*(idum-k*IQ1)-k*IR1
  if (idum.lt.0) idum=idum+IM1
  k=idum2/IQ2
  idum2=IA2*(idum2-k*IQ2)-k*IR2
  if (idum2.lt.0) idum2=idum2+IM2
  j=1+iy/NDIV
  iy=iv(j)-idum2
  iv(j)=idum
  if(iy.lt.1)iy=iy+IMM1
  ran2_dp=min(AM*iy,RNMX)
  return
  
end function ran2_dp
   
