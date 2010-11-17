!------------------------------------------------------------------------------
! ssht_inverse
!
!! Compute inverse spin spherical harmonic transform using either DH or MW 
!! sampling.  Sampled functions on the sphere and harmonic coefficients are
!! read/written to/from text files.
!!
!! Usage: ssht_inverse
!!   - [-help]: Display usage information.
!!   - [-inp filename_in]: Name of input file containing harmonic coefficients.
!!   - [-out filename_out]: Name of output file containing function samples.
!!   - [-method method]: Method (i.e. sampling) to use (DH; MW).
!!   - [-L L]: Harmonic band-limit.
!!   - [-spin spin]: Spin order. 
!!   - [-reality reality]: Integer specifying whether the function on the 
!!     sphere is real (must only be specified for spin 0 signals) 
!!     (0=false; 1=true).
!!   - [-verbosity verbosity]: Integer specify verbosity level from 0 to 5.
!     
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!
! Revisions:
!   November 2010 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program ssht_inverse

  use ssht_types_mod
  use ssht_error_mod
  use ssht_core_mod
  use s2_sky_mod

  implicit none

  character(len=*), parameter ::  METHOD_DH = 'DH'
  character(len=*), parameter ::  METHOD_MW = 'MW'
  character(len=STRING_LEN) :: method
  character(len=STRING_LEN) :: filename_in, filename_out
  character(len=*), parameter ::  IO_FORMAT = '(2e25.15)'
  integer :: L, spin
  integer :: verbosity = 0, reality = 0
  integer :: t, p, ind, ntheta
  real(dp) :: re, im
  integer :: fileid
  integer :: fail = 0
  complex(dpc), allocatable :: f(:,:), flm(:)
  real(dp), allocatable :: f_real(:,:)

  ! Parse options from command line.
  call parse_options()

  ! Check reality flag valid.
  if (spin /= 0 .and. reality == 1) then
     call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_inverse', &
          comment_add='Reality flag may only be set for spin 0 signals.')
  end if

  ! Set ntheta depending on method.
  select case(method)
     case(METHOD_DH)
        ntheta = 2*L
     case(METHOD_MW)
        ntheta = L
     case default
        call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_inverse', &
             comment_add='Invalid method.')
  end select

  ! Allocate space.
  if (reality == 1) then
     allocate(f_real(0:ntheta-1, 0:2*L-2), stat=fail)
  else
     allocate(f(0:ntheta-1, 0:2*L-2), stat=fail)
  end if
  allocate(flm(0:L**2-1), stat=fail)
  if(fail /= 0) then
     call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_inverse')
  end if  
  if (reality == 1) then
     f_real(0:ntheta-1, 0:2*L-2) = 0d0
  else
     f(0:ntheta-1, 0:2*L-2) = cmplx(0d0,0d0)
  end if
  flm(0:L**2-1) = cmplx(0d0,0d0)

  ! Read flms.
  fileid = 11
  open(fileid, file=trim(filename_in), &
       form='formatted', status='old')  
  do ind = 0,L**2 - 1
     read(fileid,IO_FORMAT) re, im
     flm(ind) = cmplx(re,im)
  end do

  ! Compute inverse transform.
  select case(method)
     case(METHOD_DH)
        if (reality == 1) then           
           call ssht_core_dh_inverse_real(f_real(0:ntheta-1, 0:2*L-2), &
                flm(0:L**2-1), L, verbosity)
        else
           call ssht_core_dh_inverse(f(0:ntheta-1, 0:2*L-2), &
                flm(0:L**2-1), L, spin, verbosity)
        end if
     case(METHOD_MW)
        if (reality == 1) then           
           call ssht_core_mw_inverse_real(f_real(0:ntheta-1, 0:2*L-2), &
                flm(0:L**2-1), L, verbosity)
        else
           call ssht_core_mw_inverse(f(0:ntheta-1, 0:2*L-2), &
                flm(0:L**2-1), L, spin, verbosity)
        end if
     case default
        call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_inverse', &
             comment_add='Invalid method.')
  end select

  ! Write f.
  fileid = 12
  open(unit=fileid, file=trim(filename_out), status='new', action='write', &
       form='formatted')
  if (reality == 1) then
     do t = 0, ntheta-1
        do p = 0, 2*L-2        
           write(fileid,IO_FORMAT) f_real(t,p)
        end do
     end do
  else
     do t = 0, ntheta-1
        do p = 0, 2*L-2        
           write(fileid,IO_FORMAT) real(f(t,p)), aimag(f(t,p))
        end do
     end do
  end if
  close(fileid)

  ! Free memory.
  if (allocated(f)) deallocate(f)
  if (allocated(f_real)) deallocate(f_real)
  deallocate(flm)

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
          write(*,'(a)') 'Usage: ssht_inverse [-inp filename_in]'
          write(*,'(a)') '                    [-out filename_out]'
          write(*,'(a)') '                    [-method method (DH; MW)]'
          write(*,'(a)') '                    [-L L]'  
          write(*,'(a)') '                    [-spin spin]'
          write(*,'(a)') '                    [-reality reality (0=false; 1=true)]'
          write(*,'(a)') '                    [-verbosity verbosity]'
          stop

       case ('-inp')
          filename_in = trim(arg)

       case ('-out')
          filename_out = trim(arg)

       case ('-method')
          method = trim(arg)

       case ('-L')
          read(arg,*) L

       case ('-spin')
          read(arg,*) spin

       case ('-reality')
          read(arg,*) reality

       case ('-verbosity')
          read(arg,*) verbosity

       case default
          print '("unknown option ",a4," ignored")', opt            

       end select
    end do

  end subroutine parse_options


end program ssht_inverse
