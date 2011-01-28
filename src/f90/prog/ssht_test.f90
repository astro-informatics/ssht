!------------------------------------------------------------------------------
! ssht_test
!
!! Applies SSHT algorithms to perform inverse and forward spherical harmonic 
!! transforms (respectively) to check that the original signal is 
!! reconstructed exactly (to numerical precision).  Test is performed on a 
!! random signal with harmonic coefficients uniformly sampled from (-1,1).
!!
!! Usage: ssht_test B spin, e.g. ssht_test 64 2
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 - November 2010
!
! Revisions:
!   November 2010 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program ssht_test

  use ssht_types_mod
  use ssht_error_mod
  use ssht_sampling_mod
  use ssht_core_mod
#ifdef NAGFOR
  use F90_UNIX_ENV
#endif

  implicit none

  interface
     subroutine ssht_test_gen_flm_real(L, flm, seed)
       use ssht_types_mod, only: dpc
       implicit none
       integer, intent(in) :: L
       complex(dpc), intent(out) :: flm(0:L**2-1)
       integer, intent(in) :: seed
     end subroutine ssht_test_gen_flm_real

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
  integer, parameter :: N_repeat = 3
  integer :: verbosity = 0
  integer :: fail = 0, seed, i_repeat
  real :: time_start, time_end

  real(dp) :: errors_dh(0:N_repeat-1)
  real :: durations_forward_dh(0:N_repeat-1)
  real :: durations_inverse_dh(0:N_repeat-1)
  real(dp) :: errors_gl(0:N_repeat-1)
  real :: durations_forward_gl(0:N_repeat-1)
  real :: durations_inverse_gl(0:N_repeat-1)
  real(dp) :: errors_mweo(0:N_repeat-1)
  real :: durations_forward_mweo(0:N_repeat-1)
  real :: durations_inverse_mweo(0:N_repeat-1)
  real(dp) :: errors_mw(0:N_repeat-1)
  real :: durations_forward_mw(0:N_repeat-1)
  real :: durations_inverse_mw(0:N_repeat-1)

  real(dp) :: errors_dh_real(0:N_repeat-1)
  real :: durations_forward_dh_real(0:N_repeat-1)
  real :: durations_inverse_dh_real(0:N_repeat-1)
  real(dp) :: errors_gl_real(0:N_repeat-1)
  real :: durations_forward_gl_real(0:N_repeat-1)
  real :: durations_inverse_gl_real(0:N_repeat-1)
  real(dp) :: errors_mweo_real(0:N_repeat-1)
  real :: durations_forward_mweo_real(0:N_repeat-1)
  real :: durations_inverse_mweo_real(0:N_repeat-1)
  real(dp) :: errors_mw_real(0:N_repeat-1)
  real :: durations_forward_mw_real(0:N_repeat-1)
  real :: durations_inverse_mw_real(0:N_repeat-1)

integer, parameter :: L = 2
  integer :: LOLD, ind, ind_check, el, el_check, m, m_check
  integer :: spin
!!$  complex(dpc), allocatable :: flm_orig(:), flm_syn(:)
!!$  complex(dpc), allocatable :: f_dh(:,:), f_gl(:,:), f_mweo(:,:), f_mw(:,:)
  real(dp), allocatable :: f_dh_real(:,:), f_gl_real(:,:), f_mweo_real(:,:), f_mw_real(:,:)
  real(dp) :: phi_sp_mweo, phi_sp_mw
  complex(dpc) :: f_sp_mweo, f_sp_mw
  real(dp) :: f_real_sp_dh, f_real_sp_gl, f_real_sp_mweo, f_real_sp_mw


  complex(dpc) :: flm_orig(0:L**2-1)
  complex(dpc) :: flm_syn(0:L**2-1)  
  complex(dpc) :: f_dh(0:2*L-1, 0:2*L-2)
  complex(dpc) :: f_gl(0:L-1, 0:2*L-2)
  complex(dpc) :: f_mweo(0:L-2, 0:2*L-2)
  complex(dpc) :: f_mw(0:L-2, 0:2*L-2)
 

  ! Initialise parameters.
  call getarg(1, arg)
  read(arg,*) LOLD
  call getarg(2, arg)
  read(arg,*) spin
  seed = 1

  ! Allocate memory.
!!$  allocate(flm_orig(0:L**2-1), stat=fail)
!!$  allocate(flm_syn(0:L**2-1), stat=fail)  
!!$  allocate(f_dh(0:2*L-1, 0:2*L-2), stat=fail)
!!$  allocate(f_gl(0:L-1, 0:2*L-2), stat=fail)
!!$  allocate(f_mweo(0:L-2, 0:2*L-2), stat=fail)
!!$  allocate(f_mw(0:L-2, 0:2*L-2), stat=fail)
  allocate(f_dh_real(0:2*L-1, 0:2*L-2), stat=fail)
  allocate(f_gl_real(0:L-1, 0:2*L-2), stat=fail)
  allocate(f_mweo_real(0:L-2, 0:2*L-2), stat=fail)
  allocate(f_mw_real(0:L-2, 0:2*L-2), stat=fail)
  if(fail /= 0) then
     call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_test')
  end if
  flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
  flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
  f_dh(0:2*L-1, 0:2*L-2) = cmplx(0d0, 0d0)
  f_gl(0:L-1, 0:2*L-2) = cmplx(0d0, 0d0)
  f_mweo(0:L-2, 0:2*L-2) = cmplx(0d0, 0d0)
  f_mw(0:L-2, 0:2*L-2) = cmplx(0d0, 0d0)
  f_dh_real(0:2*L-1, 0:2*L-2) = 0d0
  f_gl_real(0:L-1, 0:2*L-2) = 0d0
  f_mweo_real(0:L-2, 0:2*L-2) = 0d0
  f_mw_real(0:L-2, 0:2*L-2) = 0d0

  ! Write program name.
  write(*,*)
  write(*,'(a)') 'SSHT test program'
  write(*,'(a)') '==============================================================='
  write(*,*)

  ! Test index conversion between (el,m) and ind.
  ind_check = 0
  do el = 0, L-1     
     do m = -el, el       
        call ssht_sampling_elm2ind(ind, el, m)
        if (ind /= ind_check) then
           call ssht_error(SSHT_ERROR_INDEX_INVALID, 'ssht_test')
        end if
        call ssht_sampling_ind2elm(el_check, m_check, ind)
        if (el /= el_check .and. m /= m_check) then
           call ssht_error(SSHT_ERROR_INDEX_INVALID, 'ssht_test')
        end if
        ind_check = ind_check + 1
     end do
  end do
  if (ind_check /= L**2) then
     call ssht_error(SSHT_ERROR_INDEX_INVALID, 'ssht_test')
  end if

  ! Run algorithm error and timing tests.
  do i_repeat = 0,N_repeat-1

     ! If spin=0 run tests on algorithms optimised for real spin=0 signal.
!!$     if (spin == 0) then 
!!$
!!$        !=========================================================================
!!$        ! DH real spin=0
!!$        write(*,'(a,i2)') 'DH real spin=0 test no.', i_repeat
!!$        flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
!!$        flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
!!$        call ssht_test_gen_flm_real(L, flm_orig, seed)
!!$        call cpu_time(time_start)
!!$        !-------------------------------------------------------------------------
!!$        call ssht_core_dh_inverse_real(f_dh_real, flm_orig, L, verbosity)
!!$        !-------------------------------------------------------------------------
!!$        call cpu_time(time_end)
!!$        durations_inverse_dh_real(i_repeat) = time_end - time_start
!!$        call cpu_time(time_start)
!!$        !-------------------------------------------------------------------------
!!$        call ssht_core_dh_forward_real(flm_syn, f_dh_real, L, verbosity)
!!$        !-------------------------------------------------------------------------
!!$        call cpu_time(time_end)
!!$        durations_forward_dh_real(i_repeat) = time_end - time_start
!!$        errors_dh_real(i_repeat) = maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))
!!$        write(*,'(a,f40.4)') ' duration_inverse (s) =', &
!!$             durations_inverse_dh_real(i_repeat)
!!$        write(*,'(a,f40.4)') ' duration_forward (s) =', &
!!$             durations_forward_dh_real(i_repeat)
!!$        write(*,'(a,e40.5)') ' error                =', &
!!$             errors_dh_real(i_repeat)
!!$        write(*,*)
!!$
!!$        !=========================================================================
!!$        ! GL real spin=0
!!$        write(*,'(a,i2)') 'GL real spin=0 test no.', i_repeat
!!$        flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
!!$        flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
!!$        call ssht_test_gen_flm_real(L, flm_orig, seed)
!!$        call cpu_time(time_start)
!!$        !-------------------------------------------------------------------------
!!$        call ssht_core_gl_inverse_real(f_gl_real, flm_orig, L, verbosity)
!!$        !-------------------------------------------------------------------------
!!$        call cpu_time(time_end)
!!$        durations_inverse_gl_real(i_repeat) = time_end - time_start
!!$        call cpu_time(time_start)
!!$        !-------------------------------------------------------------------------
!!$        call ssht_core_gl_forward_real(flm_syn, f_gl_real, L, verbosity)
!!$        !-------------------------------------------------------------------------
!!$        call cpu_time(time_end)
!!$        durations_forward_gl_real(i_repeat) = time_end - time_start
!!$        errors_gl_real(i_repeat) = maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))
!!$        write(*,'(a,f40.4)') ' duration_inverse (s) =', &
!!$             durations_inverse_gl_real(i_repeat)
!!$        write(*,'(a,f40.4)') ' duration_forward (s) =', &
!!$             durations_forward_gl_real(i_repeat)
!!$        write(*,'(a,e40.5)') ' error                =', &
!!$             errors_gl_real(i_repeat)
!!$        write(*,*)
!!$
!!$        !=========================================================================
!!$        ! MWEO real spin=0
!!$        write(*,'(a,i2)') 'MWEO real spin=0 test no.', i_repeat
!!$        flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
!!$        flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
!!$        call ssht_test_gen_flm_real(L, flm_orig, seed)
!!$        call cpu_time(time_start)
!!$        !-------------------------------------------------------------------------
!!$        call ssht_core_mweo_inverse_real(f_mweo_real, &
!!$             f_real_sp_mweo, flm_orig, L, verbosity)
!!$        !-------------------------------------------------------------------------
!!$        call cpu_time(time_end)
!!$        durations_inverse_mweo_real(i_repeat) = time_end - time_start
!!$        call cpu_time(time_start)
!!$        !-------------------------------------------------------------------------
!!$        call ssht_core_mweo_forward_real(flm_syn, &
!!$             f_mweo_real, f_real_sp_mweo, L, verbosity)
!!$        !-------------------------------------------------------------------------
!!$        call cpu_time(time_end)
!!$        durations_forward_mweo_real(i_repeat) = time_end - time_start
!!$        errors_mweo_real(i_repeat) = maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))
!!$        write(*,'(a,f40.4)') ' duration_inverse (s) =', &
!!$             durations_inverse_mweo_real(i_repeat)
!!$        write(*,'(a,f40.4)') ' duration_forward (s) =', &
!!$             durations_forward_mweo_real(i_repeat)
!!$        write(*,'(a,e40.5)') ' error                =', &
!!$             errors_mweo_real(i_repeat)
!!$        write(*,*)
!!$
!!$        !=========================================================================
!!$        ! MW real spin=0
!!$        write(*,'(a,i2)') 'MW real spin=0 test no.', i_repeat
!!$        flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
!!$        flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
!!$        call ssht_test_gen_flm_real(L, flm_orig, seed)
!!$        call cpu_time(time_start)
!!$        !-------------------------------------------------------------------------
!!$        call ssht_core_mw_inverse_real(f_mw_real, &
!!$             f_real_sp_mw, flm_orig, L, verbosity)
!!$        !-------------------------------------------------------------------------
!!$        call cpu_time(time_end)
!!$        durations_inverse_mw_real(i_repeat) = time_end - time_start
!!$        call cpu_time(time_start)
!!$        !-------------------------------------------------------------------------
!!$        call ssht_core_mw_forward_real(flm_syn, &
!!$             f_mw_real, f_real_sp_mw, L, verbosity)
!!$        !-------------------------------------------------------------------------
!!$        call cpu_time(time_end)
!!$        durations_forward_mw_real(i_repeat) = time_end - time_start
!!$        errors_mw_real(i_repeat) = maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))
!!$        write(*,'(a,f40.4)') ' duration_inverse (s) =', &
!!$             durations_inverse_mw_real(i_repeat)
!!$        write(*,'(a,f40.4)') ' duration_forward (s) =', &
!!$             durations_forward_mw_real(i_repeat)
!!$        write(*,'(a,e40.5)') ' error                =', &
!!$             errors_mw_real(i_repeat)
!!$        write(*,*)
!!$
!!$
!!$     end if

!!$     !=========================================================================
!!$     ! DH
!!$     write(*,'(a,i2)') 'DH test no.', i_repeat
!!$     flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
!!$     flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
!!$     call ssht_test_gen_flm_complex(L, spin, flm_orig, seed)
!!$     call cpu_time(time_start)
!!$     !-------------------------------------------------------------------------
!!$     call ssht_core_dh_inverse(f_dh, flm_orig, L, spin, verbosity)
!!$     !-------------------------------------------------------------------------
!!$     call cpu_time(time_end)
!!$     durations_inverse_dh(i_repeat) = time_end - time_start
!!$     call cpu_time(time_start)
!!$     !-------------------------------------------------------------------------
!!$     call ssht_core_dh_forward(flm_syn, f_dh, L, spin, verbosity)
!!$     !-------------------------------------------------------------------------
!!$     call cpu_time(time_end)
!!$     durations_forward_dh(i_repeat) = time_end - time_start
!!$     errors_dh(i_repeat) = maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))
!!$     write(*,'(a,f40.4)') ' duration_inverse (s) =', &
!!$          durations_inverse_dh(i_repeat)
!!$     write(*,'(a,f40.4)') ' duration_forward (s) =', &
!!$          durations_forward_dh(i_repeat)
!!$     write(*,'(a,e40.5)') ' error                =', &
!!$          errors_dh(i_repeat)
!!$     write(*,*)
!!$
!!$
!!$     !=========================================================================
!!$     ! GL
!!$     write(*,'(a,i2)') 'GL test no.', i_repeat
!!$     flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
!!$     flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
!!$     call ssht_test_gen_flm_complex(L, spin, flm_orig, seed)
!!$     call cpu_time(time_start)
!!$     !-------------------------------------------------------------------------
!!$     call ssht_core_gl_inverse(f_gl, flm_orig, L, spin, verbosity)
!!$     !-------------------------------------------------------------------------
!!$     call cpu_time(time_end)
!!$     durations_inverse_gl(i_repeat) = time_end - time_start
!!$     call cpu_time(time_start)
!!$     !-------------------------------------------------------------------------
!!$     call ssht_core_gl_forward(flm_syn, f_gl, L, spin, verbosity)
!!$     !-------------------------------------------------------------------------
!!$     call cpu_time(time_end)
!!$     durations_forward_gl(i_repeat) = time_end - time_start
!!$     errors_gl(i_repeat) = maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))
!!$     write(*,'(a,f40.4)') ' duration_inverse (s) =', &
!!$          durations_inverse_gl(i_repeat)
!!$     write(*,'(a,f40.4)') ' duration_forward (s) =', &
!!$          durations_forward_gl(i_repeat)
!!$     write(*,'(a,e40.5)') ' error                =', &
!!$          errors_gl(i_repeat)
!!$     write(*,*)
!!$
!!$     !=========================================================================
!!$     ! MWEO
!!$     write(*,'(a,i2)') 'MWEO test no.', i_repeat
!!$     flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
!!$     flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
!!$     call ssht_test_gen_flm_complex(L, spin, flm_orig, seed)
!!$     call cpu_time(time_start)
!!$     !-------------------------------------------------------------------------
!!$     call ssht_core_mweo_inverse(f_mweo, f_sp_mweo, phi_sp_mweo, &
!!$          flm_orig, L, spin, verbosity)
!!$     !-------------------------------------------------------------------------
!!$     call cpu_time(time_end)
!!$     durations_inverse_mweo(i_repeat) = time_end - time_start
!!$     call cpu_time(time_start)
!!$     !-------------------------------------------------------------------------
!!$     call ssht_core_mweo_forward(flm_syn, f_mweo, &
!!$          f_sp_mweo, phi_sp_mweo, L, spin, verbosity)
!!$     !-------------------------------------------------------------------------
!!$     call cpu_time(time_end)
!!$     durations_forward_mweo(i_repeat) = time_end - time_start
!!$     errors_mweo(i_repeat) = maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))
!!$     write(*,'(a,f40.4)') ' duration_inverse (s) =', &
!!$          durations_inverse_mweo(i_repeat)
!!$     write(*,'(a,f40.4)') ' duration_forward (s) =', &
!!$          durations_forward_mweo(i_repeat)
!!$     write(*,'(a,e40.5)') ' error                =', &
!!$          errors_mweo(i_repeat)
!!$     write(*,*)

     !=========================================================================
     ! MW
     write(*,'(a,i2)') 'MW test no.', i_repeat
     flm_orig(0:L**2-1) = cmplx(0d0, 0d0)
     flm_syn(0:L**2-1) = cmplx(0d0, 0d0)
     call ssht_test_gen_flm_complex(L, spin, flm_orig, seed)
do el = 0,L*L-1
   write(*,*) 'flm_orig[el]=', flm_orig(el)
end do

     call cpu_time(time_start)
     !-------------------------------------------------------------------------
     call ssht_core_mw_inverse(f_mw, f_sp_mw, phi_sp_mw, &
          flm_orig, L, spin, verbosity)
     !-------------------------------------------------------------------------
     call cpu_time(time_end)
     durations_inverse_mw(i_repeat) = time_end - time_start
     call cpu_time(time_start)
     !-------------------------------------------------------------------------
     call ssht_core_mw_forward(flm_syn, f_mw, &
          f_sp_mw, phi_sp_mw, L, spin, verbosity)
     !-------------------------------------------------------------------------
     call cpu_time(time_end)
     durations_forward_mw(i_repeat) = time_end - time_start
     errors_mw(i_repeat) = maxval(abs(flm_orig(0:L**2-1) - flm_syn(0:L**2-1)))
     write(*,'(a,f40.4)') ' duration_inverse (s) =', &
          durations_inverse_mw(i_repeat)
     write(*,'(a,f40.4)') ' duration_forward (s) =', &
          durations_forward_mw(i_repeat)
     write(*,'(a,e40.5)') ' error                =', &
          errors_mw(i_repeat)
     write(*,*)

  end do

  ! Print summary.
  write(*,'(a)') '==============================================================='
  write(*,'(a)') 'Summary'
  write(*,*)
  write(*,'(a,i40)') 'N_repeat              =', N_repeat
  write(*,'(a,i40)') 'L                     =', L
  write(*,'(a,i40)') 'spin                  =', spin
  write(*,*)

  if (spin == 0) then
     write(*,'(a,i2)') 'DH real'
     write(*,'(a,f30.5)') ' Average forward transform time =', &
          sum(durations_forward_dh_real(0:N_repeat-1)) / real(N_repeat)
     write(*,'(a,f30.5)') ' Average inverse transform time =', &
          sum(durations_inverse_dh_real(0:N_repeat-1)) / real(N_repeat)
     write(*,'(a,e30.5)') ' Average max error              =', &
          sum(errors_dh_real(0:N_repeat-1)) / real(N_repeat)
     write(*,*)

     write(*,'(a,i2)') 'GL real'
     write(*,'(a,f30.5)') ' Average forward transform time =', &
          sum(durations_forward_gl_real(0:N_repeat-1)) / real(N_repeat)
     write(*,'(a,f30.5)') ' Average inverse transform time =', &
          sum(durations_inverse_gl_real(0:N_repeat-1)) / real(N_repeat)
     write(*,'(a,e30.5)') ' Average max error              =', &
          sum(errors_gl_real(0:N_repeat-1)) / real(N_repeat)
     write(*,*)

     write(*,'(a,i2)') 'MWEO real'
     write(*,'(a,f30.5)') ' Average forward transform time =', &
          sum(durations_forward_mweo_real(0:N_repeat-1)) / real(N_repeat)
     write(*,'(a,f30.5)') ' Average inverse transform time =', &
          sum(durations_inverse_mweo_real(0:N_repeat-1)) / real(N_repeat)
     write(*,'(a,e30.5)') ' Average max error              =', &
          sum(errors_mweo_real(0:N_repeat-1)) / real(N_repeat)
     write(*,*)

     write(*,'(a,i2)') 'MW real'
     write(*,'(a,f30.5)') ' Average forward transform time =', &
          sum(durations_forward_mw_real(0:N_repeat-1)) / real(N_repeat)
     write(*,'(a,f30.5)') ' Average inverse transform time =', &
          sum(durations_inverse_mw_real(0:N_repeat-1)) / real(N_repeat)
     write(*,'(a,e30.5)') ' Average max error              =', &
          sum(errors_mw_real(0:N_repeat-1)) / real(N_repeat)
     write(*,*)
  end if

  write(*,'(a,i2)') 'DH'
  write(*,'(a,f30.5)') ' Average forward transform time =', &
       sum(durations_forward_dh(0:N_repeat-1)) / real(N_repeat)
  write(*,'(a,f30.5)') ' Average inverse transform time =', &
       sum(durations_inverse_dh(0:N_repeat-1)) / real(N_repeat)
  write(*,'(a,e30.5)') ' Average max error              =', &
       sum(errors_dh(0:N_repeat-1)) / real(N_repeat)
  write(*,*)

  write(*,'(a,i2)') 'GL'
  write(*,'(a,f30.5)') ' Average forward transform time =', &
       sum(durations_forward_gl(0:N_repeat-1)) / real(N_repeat)
  write(*,'(a,f30.5)') ' Average inverse transform time =', &
       sum(durations_inverse_gl(0:N_repeat-1)) / real(N_repeat)
  write(*,'(a,e30.5)') ' Average max error              =', &
       sum(errors_gl(0:N_repeat-1)) / real(N_repeat)
  write(*,*)

  write(*,'(a,i2)') 'MWEO'
  write(*,'(a,f30.5)') ' Average forward transform time =', &
       sum(durations_forward_mweo(0:N_repeat-1)) / real(N_repeat)
  write(*,'(a,f30.5)') ' Average inverse transform time =', &
       sum(durations_inverse_mweo(0:N_repeat-1)) / real(N_repeat)
  write(*,'(a,e30.5)') ' Average max error              =', &
       sum(errors_mweo(0:N_repeat-1)) / real(N_repeat)
  write(*,*)

  write(*,'(a,i2)') 'MW'
  write(*,'(a,f30.5)') ' Average forward transform time =', &
       sum(durations_forward_mw(0:N_repeat-1)) / real(N_repeat)
  write(*,'(a,f30.5)') ' Average inverse transform time =', &
       sum(durations_inverse_mw(0:N_repeat-1)) / real(N_repeat)
  write(*,'(a,e30.5)') ' Average max error              =', &
       sum(errors_mw(0:N_repeat-1)) / real(N_repeat)
  write(*,*)

  ! Deallocate memory.
!!$  deallocate(flm_orig, flm_syn)
!!$  deallocate(f_dh, f_gl, f_mweo, f_mw)
!!$  deallocate(f_dh_real, f_gl_real, f_mweo_real, f_mw_real)

end program ssht_test


!--------------------------------------------------------------------------
! ssht_test_gen_flm_real
!
!! Generate random spherical harmonic coefficients of a real spin=0
!! signal.
!!
!! Variables:
!!   - L: Harmonic band-limit [input].
!!   - flm(0:L**2-1): Random spherical harmonic coefficients generated 
!!     [output].
!!   - seed: Integer seed required for random number generator [input].
!
!! @author J. D. McEwen
!! @version 0.1 October 2010
! 
! Revisions:
!   October 2010 - Jason McEwen
!--------------------------------------------------------------------------

subroutine ssht_test_gen_flm_real(L, flm, seed)

  use ssht_types_mod, only: dpc
  use ssht_sampling_mod

  implicit none

  interface 
     function ran2_dp(idum)
       use ssht_types_mod, only: dp
       real(dp) :: ran2_dp
       integer :: idum
     end function ran2_dp
  end interface

  integer, intent(in) :: L
  complex(dpc), intent(out) :: flm(0:L**2-1)
  integer, intent(in) :: seed

  integer :: el, m, ind
  complex(dpc) :: tmp

  flm(0:L**2-1) = cmplx(0d0, 0d0)
  do el = 0,L-1
     m = 0
     call ssht_sampling_elm2ind(ind, el, m)  
     tmp = cmplx(2d0*ran2_dp(seed)-1d0, 0d0)
     flm(ind) = tmp

     do m = 1,el
        call ssht_sampling_elm2ind(ind, el, m)  
        tmp = cmplx(2d0*ran2_dp(seed)-1d0, 2d0*ran2_dp(seed)-1d0)
        flm(ind) = tmp
        call ssht_sampling_elm2ind(ind, el, -m)  
        flm(ind) = (-1)**m * conjg(tmp)
     end do

  end do

end subroutine ssht_test_gen_flm_real


!--------------------------------------------------------------------------
! ssht_test_gen_flm_complex
!
!! Generate random spherical harmonic coefficients of a complex signal.
!!
!! Variables:
!!   - L: Harmonic band-limit [input].
!!   - flm(0:L**2-1): Random spherical harmonic coefficients generated 
!!     [output].
!!   - seed: Integer seed required for random number generator [input].
!
!! @author J. D. McEwen
!! @version 0.1 October 2010
! 
! Revisions:
!   October 2010 - Jason McEwen
!--------------------------------------------------------------------------

subroutine ssht_test_gen_flm_complex(L, spin, flm, seed)

  use ssht_types_mod, only: dpc
  use ssht_sampling_mod

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

  flm(0:L**2-1) = cmplx(0d0, 0d0)

  call ssht_sampling_elm2ind(ind_lo, abs(spin), 0)
  do ind = ind_lo,L**2 - 1
!     call ssht_core_ind2elm(el, m, ind)  
!     if(el >= abs(spin)) then
        flm(ind) = cmplx(2d0*ran2_dp(seed)-1d0, 2d0*ran2_dp(seed)-1d0)
!     end if
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

