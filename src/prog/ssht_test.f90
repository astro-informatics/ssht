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

     subroutine ssht_test_gen_flm_complex(L, spin, flm, seed)
       use ssht_types_mod, only: dpc
       implicit none
       integer, intent(in) :: L
       integer, intent(in) :: spin
       complex(dpc), intent(out) :: flm(0:L**2-1)
       integer, intent(in) :: seed
     end subroutine ssht_test_gen_flm_complex

  end interface

  character(len=64) :: arg
  integer, parameter :: N_repeat = 1
  integer :: fail = 0, seed, i_repeat
  real(dp) :: error_flm(0:N_repeat-1)
  logical :: admiss_pass
  real :: time_start, time_end
  real :: durations_analysis(0:N_repeat-1)
  real :: durations_synthesis(0:N_repeat-1)

  integer :: L, ind, ind_check, el, el_check, m, m_check
  integer :: spin
  complex(dpc), allocatable :: flm_orig(:), flm_syn(:)
  complex(dpc), allocatable :: f_dh(:,:), f_hw(:,:), f_mweo(:,:)

  ! Initialise parameters.
  call getarg(1, arg)
  read(arg,*) L
  call getarg(2, arg)
  read(arg,*) spin
  seed = 1

  ! Allocate memory.
  allocate(flm_orig(0:L**2-1), stat=fail)
  allocate(flm_syn(0:L**2-1), stat=fail)  
  allocate(f_dh(0:2*L-1, 0:2*L-2), stat=fail)
  allocate(f_hw(0:L-1, 0:2*L-1), stat=fail)
  allocate(f_mweo(0:L-1, 0:2*L-2), stat=fail)
  if(fail /= 0) then
     call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_test')
  end if
  flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
  flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
  f_dh(0:2*L-1, 0:2*L-2) = cmplx(0d0, 0d0)
  f_hw(0:L-1, 0:2*L-1) = cmplx(0d0, 0d0)
  f_mweo(0:L-1, 0:2*L-2) = cmplx(0d0, 0d0)



  write(*,*)
  write(*,'(a)') 'SSHT test program'
  write(*,'(a)') '================='
  write(*,*)




  ! Test conversion between (el,m) and ind.
  !L = B
  ind_check = 0
  do el = 0, L-1     
     do m = -el, el       
!write(*,*) 'Testing el,m'
        call ssht_core_elm2ind(ind, el, m)
        if (ind /= ind_check) then
stop "failed"
        end if
        call ssht_core_ind2elm(el_check, m_check, ind)
        if (el /= el_check .and. m /= m_check) then
stop "failed"
        end if
        ind_check = ind_check + 1
     end do
  end do

write(*,*) 'ind_check = ', ind_check
write(*,*) 'L**2 = ', L**2



  do i_repeat = 0,N_repeat-1

     ! Generate harmonic coefficients of random test signal.
!     call ssht_test_gen_flm(B-1, flm_orig, seed)


     flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
     flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
     call ssht_test_gen_flm_complex(L, spin, flm_orig, seed)
     !call ssht_core_dh_inverse_direct(f_dh, flm2_orig, L, spin)
     call ssht_core_dh_inverse_sov(f_dh, flm_orig, L, spin)
     call ssht_core_dh_forward_sov(flm_syn, f_dh, L, spin)

     write(*,'(a,e43.5)') 'HERE IT IS, MAXERR: ', maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))


!!$do ind = 0,L*2-1
!!$   write(*,'(a,2f10.5)') 'flm_orig ', flm_orig(ind)
!!$end do
!!$write(*,*)
!!$do ind = 0,L*2-1
!!$   write(*,'(a,2f10.5)') 'flm_syn ', flm_syn(ind)
!!$end do

     flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
     flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
     call ssht_test_gen_flm_complex(L, spin, flm_orig, seed)
     call ssht_core_hw_inverse_direct(f_hw, flm_orig, L, spin)
     call ssht_core_hw_forward_direct(flm_syn, f_hw, L, spin)

     write(*,'(a,e43.5)') 'HERE IT IS, MAXERR: ', maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))

     flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
     flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
     call ssht_test_gen_flm_complex(L, spin, flm_orig(0:L**2-1), seed)
     !call ssht_core_mweo_inverse_direct(f_mweo, flm2_orig, L, spin)
     !call ssht_core_mweo_inverse_sov_direct(f_mweo, flm2_orig, L, spin)
     call ssht_core_mweo_inverse_sov(f_mweo(0:L-1, 0:2*L-2), flm_orig(0:L**2-1), L, spin)
     call ssht_core_mweo_forward_sov_conv(flm_syn, f_mweo, L, spin)
     !call ssht_core_mw_forward_direct(flm_syn, f_mweo, L, spin)

     write(*,'(a,e43.5)') 'HERE IT IS, MAXERR: ', maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))

     flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
     flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
     call ssht_test_gen_flm_complex(L, spin, flm_orig(0:L**2-1), seed)
     !call ssht_core_mweo_inverse_direct(f_mweo, flm2_orig, L, spin)
     !call ssht_core_mweo_inverse_sov_direct(f_mweo, flm2_orig, L, spin)
     call ssht_core_mw_inverse_sov_direct(f_mweo(0:L-1, 0:2*L-2), flm_orig(0:L**2-1), L, spin)
     !call ssht_core_mweo_forward_sov_conv(flm_syn, f_mweo, L, spin)
     call ssht_core_mw_forward_direct(flm_syn, f_mweo, L, spin)

     write(*,'(a,e43.5)') 'HERE IT IS, MAXERR: ', maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))



!deallocate(flm2_orig, flm2_syn, f_dh, f_mw, f_mweo)


     ! Compute kernels, scaling function and directionality coefficients.
     write(*,'(a,i2)') 'Initialisation no.', i_repeat
     call cpu_time(time_start)
     !...
     call cpu_time(time_end)
     write(*,'(a,f43.2)') ' duration =', time_end - time_start


     ! Compute wavelet and scaling coefficients.
     write(*,'(a,i2)') 'Analysis no.', i_repeat
     call cpu_time(time_start)
     !...     
     call cpu_time(time_end)
     durations_analysis(i_repeat) = time_end - time_start
     write(*,'(a,f43.2)') ' duration =', durations_analysis(i_repeat)

     ! Synthesis harmonic coefficients of signal from wavelet and scaling 
     ! coefficients.
     write(*,'(a,i2)') 'Synthesis no.', i_repeat
     call cpu_time(time_start)
     !...     
     call cpu_time(time_end)
     durations_synthesis(i_repeat) = time_end - time_start
     write(*,'(a,f43.2)') ' duration =', durations_synthesis(i_repeat)

     !error_flm(i_repeat) = maxval(abs(flm_orig(0:B-1,0:B-1) - flm_syn(0:B-1,0:B-1)))
     write(*,'(a,e43.5)') ' error =   ', error_flm(i_repeat)
     write(*,*)

  end do

!!$  write(*,'(a)') 'Summary:'
!!$  write(*,'(a,i30)') 'N_repeat               =', N_repeat
!!$  write(*,'(a,i30)') 'B                      =', B
!!$  write(*,'(a,i30)') 'N                      =', N
!!$  write(*,'(a,i30)') 'J                      =', J	
!!$  write(*,'(a,f30.5)') 'alpha                  =', alpha
!!$  write(*,'(a,f30.5)') 'Average analysis time  =', sum(durations_analysis(0:N_repeat-1)) / real(N_repeat)
!!$  write(*,'(a,f30.5)') 'Average synthesis time =', sum(durations_synthesis(0:N_repeat-1)) / real(N_repeat)
!!$  write(*,'(a,e30.5)') 'Average max error      =', sum(error_flm(0:N_repeat-1)) / real(N_repeat)
!!$  write(*,*)
!!$
!!$  if( abs(sum(error_flm(0:N_repeat-1)) / real(N_repeat)) < 1d-6) then
!!$     write(*,'(a)') 'Tests passed!'
!!$  else
!!$     write(*,'(a)') 'Tests failed!'
!!$  end if
!!$  write(*,*) 

  ! Deallocate memory.
  deallocate(flm_orig, flm_syn)
  deallocate(f_dh, f_hw, f_mweo)

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


subroutine ssht_test_gen_flm_complex(L, spin, flm, seed)

  use ssht_types_mod, only: dpc
  use ssht_core_mod

  implicit none

  interface 
     function ran2_dp(idum)
       use ssht_types_mod, only: dp
       real(dp) :: ran2_dp
       integer :: idum
     end function ran2_dp
  end interface

  integer, intent(in) :: L
  integer, intent(in) :: spin
  complex(dpc), intent(out) :: flm(0:L**2-1)
  integer, intent(in) :: seed

  integer :: ind, ind_lo, el, m

  flm(0:L**2-1) = 0d0

  !call ssht_core_elm2ind(ind_lo, abs(spin), 0)
  do ind = 0,L**2 - 1
     call ssht_core_ind2elm(el, m, ind)  
     if(el >= abs(spin)) then
        flm(ind) = cmplx(2d0*ran2_dp(seed)-1d0, 2d0*ran2_dp(seed)-1d0)
     end if
  end do

end subroutine ssht_test_gen_flm_complex


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

