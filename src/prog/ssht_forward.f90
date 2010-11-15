
program ssht_forward

  use ssht_types_mod
  use ssht_error_mod
  use ssht_fileio_mod
  use ssht_core_mod
  use s2_sky_mod

  implicit none

  character(len=*), parameter ::  METHOD_DH = 'DH'
  character(len=*), parameter ::  METHOD_MW = 'MW'
  character(len=STRING_LEN) :: method
  character(len=STRING_LEN) :: filename_in, filename_out
  character(len=*), parameter ::  IO_FORMAT = '(2e25.15)'
  integer :: L, spin
  integer :: t, p, ind, ntheta
  real(dp) :: re, im
  integer :: fileid
  integer :: fail = 0
  complex(dpc), allocatable :: f(:,:), flm(:)

  ! Parse options from command line.
  call parse_options()

  ! Set ntheta depending on method.
  select case(method)
     case(METHOD_DH)
        ntheta = 2*L
     case(METHOD_MW)
        ntheta = L
     case default
        call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_forward', &
             comment_add='Invalid method.')
  end select

  ! Allocate space.
  allocate(f(0:ntheta-1, 0:2*L-2), stat=fail)
  allocate(flm(0:L**2-1), stat=fail)
  if(fail /= 0) then
     call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_forward')
  end if  
  f(0:ntheta-1, 0:2*L-2) = cmplx(0d0,0d0)
  flm(0:L**2-1) = cmplx(0d0,0d0)

  ! Read function.
  fileid = 11
  open(fileid, file=trim(filename_in), &
       form='formatted', status='old')  
  do t = 0, ntheta-1
     do p = 0, 2*L-2   
        read(fileid,IO_FORMAT) re, im
        f(t,p) = cmplx(re,im)
     end do
  end do
  close(fileid)

  ! Compute forward transform.
  select case(method)
     case(METHOD_DH)
        call ssht_core_dh_forward_sov(flm(0:L**2-1), &
             f(0:ntheta-1, 0:2*L-2), L, spin)
     case(METHOD_MW)
        call ssht_core_mweo_forward_sov_conv(flm(0:L**2-1), &
             f(0:ntheta-1, 0:2*L-2), L, spin)

!        call ssht_core_mw_forward_direct(flm(0:L**2-1), &
!             f(0:ntheta-1, 0:2*L-2), L, spin)

     case default
        call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_forward', &
             comment_add='Invalid method.')
  end select

  ! Write flms.
  fileid = 12
  open(unit=fileid, file=trim(filename_out), status='new', action='write', &
       form='formatted')
  do ind = 0,L**2 - 1
     write(fileid,IO_FORMAT) real(flm(ind)), aimag(flm(ind))
  end do
  close(fileid)

  ! Free memory.
  deallocate(f, flm)

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
          write(*,'(a)') 'Usage: ssht_forward [-inp filename_in]'
          write(*,'(a)') '                    [-out filename_out]'
          write(*,'(a)') '                    [-method method (DH; MW)]'
          write(*,'(a)') '                    [-L L]'  
          write(*,'(a)') '                    [-spin spin]'
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

       case default
          print '("unknown option ",a4," ignored")', opt            

       end select
    end do

  end subroutine parse_options


end program ssht_forward
