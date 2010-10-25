!------------------------------------------------------------------------------
! ssht_dl_mod
!
!! Functionality to compute specified plane of the Wigner dl matrix.
!!
!! Note:
!!  - Copied from s2 library so SSHT library is self contained.
!
!! @author D. J. Mortlock
!! @version 0.1 Compiled in SSHT library by JDM October 2007
!
! Revisions:
!   November 2007 - Jason McEwen
!------------------------------------------------------------------------------

module ssht_dl_mod

  use ssht_types_mod, only: dp

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    ssht_dl_beta_operator


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! ssht_dl_beta_operator
    !
    !! Calculates the lth plane of a d-matrix using Turok & Bucher's
    !! operator-based method.
    !!
    !! Variables:
    !!   - dl: Dl matrix values computed for lth plane.
    !!   - beta: Beta euler angle to compute dl matrix for.
    !!   - l: Plane of dl matrix to compute values for.
    !
    !! @author D. J. Mortlock
    !! @version 0.1 Compiled in ssht library by JDM October 2007
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock 
    !--------------------------------------------------------------------------
    
    subroutine ssht_dl_beta_operator(dl, beta, l)

      integer, intent(in) :: l
      real(dp), intent(out) :: dl(-l:l,-l:l)
      real(dp), intent(in) :: beta
  
      ! Fill in the top quarter of the d-matrix.
      call ssht_dl_beta_operator_core(dl(-l:l,-l:l), beta, l)
  
      ! Use its symmetry properties to fill in the rest of it.
      call ssht_dl_beta_operator_fill(dl(-l:l,-l:l), l)
  
    end subroutine ssht_dl_beta_operator


    !--------------------------------------------------------------------------
    ! ssht_dl_beta_operator_core
    !
    !! Does the left quarter of the d-matrix. Beta is the angle of rotation,
    !! l is the plane of the matrix required and dl is the two dimensional 
    !! array representing the lth plane of the matrix. This array must be
    !! already allocated on entry with dimensions given by the command
    !! allocate(dl(- lmax: lmax, - lmax: lmax)) where lmax >= l.
    !!
    !! Variables:
    !!   - dl: Dl matrix values computed for lth plane.
    !!   - beta: Beta euler angle to compute dl matrix for.
    !!   - l: Plane of dl matrix to compute values for.
    !
    !! @author D. J. Mortlock
    !! @version 0.1 Compiled in ssht library by JDM September 2005
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock 
    !--------------------------------------------------------------------------

    subroutine ssht_dl_beta_operator_core(dl, beta, l)
      integer, intent(in) :: l
      real(kind = dp), intent(out) :: dl(-l:l,-l:l)
      real(kind = dp), intent(in) :: beta

  
      real(kind = dp) :: lambda, xllp1, xm, big, bigi, bigi2, c, s, omc, &
        lbig, expo, renorm, ratio, c2, t, si, ang, big_const
      real(kind = dp), pointer :: cp(:) => null(), cpi(:) => null(), &
        cp2(:) => null(), log_first_row(:) => null(), sign(:) => null(), &
        lrenorm(:) => null()
      integer :: index, i, m, im, j, lp1
  
      ! Added by Jason McEwen 19/11/05
      ! If beta=0 then dl=identity
      real(kind = dp) :: ZERO_TOL = 1d-5
      if(abs(beta) < ZERO_TOL) then
        dl(-l:l,-l:l) = 0d0
        do i = -l,l
          dl(i,i) = 1d0
        end do
        return      ! ** Exit routine
      end if
  
      lp1 = l + 1
      big_const = 1.0d150
  
      allocate(cp(2 * l + 1))
      allocate(cpi(2 * l + 1))
      allocate(cp2(2 * l + 1))
      allocate(log_first_row(2 * l + 1))
      allocate(sign(2 * l + 1))
      allocate(lrenorm(2 * l + 1))

      ang=beta
      c=dcos(ang)
      s=dsin(ang)
      si=1.d0/s
      t=dtan(-ang/2.d0)
      c2=dcos(ang/2.d0)
      omc=1.d0-c
  
      do i=1,2*l+1
        lrenorm(i)=0.d0
      end do
  
      ! Compute first row.
  
      log_first_row(1)=(2.d0*real(l,kind=dp))*dlog(dabs(c2))
      sign(1)=1.d0
      do i=2, 2*l+1
        m=l+1-i
        ratio=dsqrt(real(l+m+1,kind=dp)/real(l-m,kind=dp))
        log_first_row(i)=log_first_row(i-1) &
          +dlog(ratio)+dlog(dabs(t))
        sign(i)=sign(i-1)*t/dabs(t)
      end do
  
      big=big_const
      lbig=dlog(big)
      bigi=1.d0/big_const
      bigi2=1.d0/big_const**2
      xllp1=real(l*(l+1),kind=dp)
  
      ! Initialising coefficients cp(m)= cplus(l-m).
  
      do m=1,l
        xm=real(l-m,kind=dp)
        cpi(m)=2.d0/dsqrt(xllp1-xm*(xm+1))
        cp(m)=1.d0/cpi(m)
      end do
      do m=2,l
        cp2(m)=cpi(m)*cp(m-1)
      end do
      dl(1 - lp1, 1 - lp1)=1.d0
      dl(2*l+1 - lp1, 1 - lp1)=1.d0
  
      ! Use recurrrence relation to fill columns down to diagonals.
  
      do index= 2,l+1
        dl(index - lp1, 1 - lp1)=1.d0
        lambda=(real(l+1,kind=dp)*omc-real(index,kind=dp)+c)*si
        dl(index - lp1, 2 - lp1)=lambda*dl(index - lp1, 1 - lp1)*cpi(1)
        if(index.gt.2) then
          do m=2,index-1
            lambda=(real(l+1,kind=dp)*omc &
              -real(index,kind=dp)+real(m,kind=dp)*c)/s
            dl(index - lp1, m+1 - lp1)= &
              lambda*cpi(m)*dl(index - lp1, m - lp1)-cp2(m) &
              *dl(index - lp1, m-1 - lp1)
            if(dl(index - lp1, m+1 - lp1).gt.big) then
              lrenorm(index)=lrenorm(index)-lbig
              do im=1,m+1
                dl(index - lp1, im - lp1)=dl(index - lp1, im - lp1)*bigi
              end do
            end if
          end do
        end if
      end do   
  
      ! Other half of triangle.
  
      do index= l+2,2*l
        dl(index - lp1, 1 - lp1)=1.d0
        lambda=(real(l+1,kind=dp)*omc-real(index,kind=dp)+c)/s
        dl(index - lp1, 2 - lp1)=lambda*dl(index - lp1, 1 - lp1)*cpi(1)
        if(index.lt.2*l) then
          do m=2,2*l-index+1
            lambda=(real(l+1,kind=dp)*omc-real(index,kind=dp) &
              +real(m,kind=dp)*c)/s
            dl(index - lp1, m+1 - lp1)= &
              lambda*cpi(m)*dl(index - lp1, m - lp1)-cp2(m)&
              *dl(index - lp1, m-1 - lp1)
            if(dl(index - lp1, m+1 - lp1).gt.big) then
              do im=1,m+1
                dl(index - lp1, im - lp1)=dl(index - lp1, im - lp1)*bigi
              end do
              lrenorm(index)=lrenorm(index)-lbig
            end if
          end do
        end if
      end do   
  
      do i=1, l+1
        renorm=sign(i)*dexp(log_first_row(i)-lrenorm(i))
        do j=1, i
          dl(i - lp1, j - lp1)= dl(i - lp1, j - lp1)*renorm
        end do
      end do
      do i=l+2,2*l+1
        expo=log_first_row(i)-lrenorm(i)
        renorm=sign(i)*dexp(log_first_row(i)-lrenorm(i))
        do j=1,2*l+2-i
          dl(i - lp1, j - lp1)=dl(i - lp1, j - lp1)*renorm
        end do
      end do
  
      ! Free memory.
      ! (Added by JDM 07/06/2004)
      deallocate(cp)
      deallocate(cpi)
      deallocate(cp2)
      deallocate(log_first_row)
      deallocate(sign)
      deallocate(lrenorm)
  
    end subroutine ssht_dl_beta_operator_core

    
    !--------------------------------------------------------------------------
    ! ssht_dl_beta_operator_fill
    !
    !! Computes the three remaining (top, bottom and right) quarters of the
    !! d-matrix, given the left quarter.
    !!
    !! Variables:
    !!   - dl: Dl matrix values computed for lth plane.
    !!   - beta: Beta euler angle to compute dl matrix for.
    !
    !! @author D. J. Mortlock
    !! @version 0.1 Compiled in ssht library by JDM September 2005
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock 
    !--------------------------------------------------------------------------

    subroutine ssht_dl_beta_operator_fill(dl, l)

      integer, intent(in) :: l
      real(kind = dp), intent(inout) :: dl(-l:l,-l:l)
  
      integer :: i, j,sgn, lp1
  
      lp1 = l + 1
  
      ! Reflect across anti-diagonal.
  
      do i = 1, l
        do j = l + 1, 2 * l + 1 - i
          dl(2 * l + 2 - i - lp1, 2 * l + 2 - j - lp1) = dl(j - lp1, i - lp1)
        end do
      end do
  
      ! Reflect across diagonal.
  
      do i = 1, l + 1
        sgn = - 1
        do j = i + 1, l + 1
          dl(i - lp1, j - lp1) = dl(j - lp1, i - lp1) * sgn
          sgn = sgn * (- 1)
        end do
      end do
  
      ! Fill in right quarter of matrix.
  
      do i = l + 2, 2 * l + 1
        sgn = (- 1)**(i + 1)
        do j = 1, 2 * l + 2 - i
          dl(j - lp1, i - lp1) = dl(i - lp1, j - lp1) * sgn
          sgn = sgn * (- 1)
        end do
        do j = i, 2 * l + 1
          dl(j - lp1, i - lp1) = dl(2 * l + 2 - i - lp1, 2 * l + 2 - j - lp1)
        end do
      end do
  
      do i = l + 2, 2 * l + 1
        do j = 2 * l + 3 - i, i - 1
          dl(j - lp1, i - lp1) = dl(2 * l + 2 - i - lp1, 2 * l + 2 - j - lp1)
        end do
      end do
  
    end subroutine ssht_dl_beta_operator_fill


end module ssht_dl_mod
