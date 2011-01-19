
program ssht_trapani

  use ssht_types_mod
  use ssht_error_mod
  use ssht_sampling_mod
  use ssht_core_mod
  !use F90_UNIX_ENV





  use ssht_dl_mod




  implicit none


  character(len=64) :: arg
  integer :: fail = 0
  real :: time_start, time_end
  integer :: el, m, mm, L
  real(dp), allocatable :: dl(:,:)
  real(dp), allocatable :: sqrt_tbl(:)

  ! Initialise parameters.
  call getarg(1, arg)
  read(arg,*) L

  write(*,*) 'Compute dl plane for L = ', L, '...'


  ! Allocate memory.
  allocate(dl(0:L-1, 0:L-1), stat=fail)
  allocate(sqrt_tbl(0:2*L+1), stat=fail)
  if(fail /= 0) then
     call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_test')
  end if
  dl(0:L-1, 0:L-1) = cmplx(0d0, 0d0)
  sqrt_tbl(0:2*(L-1)+1) = 0d0

  ! Compute square root table.
  do el = 0, 2*(L-1) + 1
     sqrt_tbl(el) = dsqrt(real(el,kind=dp))
  end do

  ! Compute quarter of dl plane recursively up to L-1
  ! (compute eighth and fill to quarter using symmetry).
  call cpu_time(time_start)
  do el = 0, L-1
     call ssht_dl_halfpi_trapani_eighth_table(dl(0:el,0:el), el, sqrt_tbl(0:2*el+1))
     call ssht_dl_halfpi_trapani_fill_eighth2quarter(dl(0:el,0:el), el)
  end do
  call cpu_time(time_end)

  write(*,*) 'Duration = ', time_end - time_start, ' seconds'

  ! Free memory.
  deallocate(dl)
  deallocate(sqrt_tbl)

end program ssht_trapani
