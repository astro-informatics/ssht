!------------------------------------------------------------------------------
! ssht_dl_mod
!
!! Functionality to compute specified plane of the Wigner dl matrix.
!!
!! Note:
!!  - Copied from s2 library so SSHT library is self contained.
!
!! @author D. J. Mortlock
!
! Revisions:
!   Ortober 2010 - Jason McEwen
!------------------------------------------------------------------------------

module ssht_dl_mod

  use ssht_types_mod, only: dp, SQRT2, PION2
  use ssht_error_mod

  implicit none

  private

  
  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
       ssht_dl_beta_operator, &
       ssht_dl_beta_risbo_full, &
       ssht_dl_beta_risbo_full_table, &
       ssht_dl_beta_risbo_half_table, &
       ssht_dl_beta_risbo_fill_half2full, &
       ssht_dl_halfpi_risbo_eighth_table, &
       ssht_dl_halfpi_risbo_fill_eighth2full, &
       ssht_dl_beta_risbo_healpix_half_table, &
       ssht_dl_halfpi_trapani_eighth, &
       ssht_dl_halfpi_trapani_eighth_table, &
       ssht_dl_halfpi_trapani_fill_eighth2full, &
       ssht_dl_halfpi_trapani_fill_eighth2quarter


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! Operator method
    !--------------------------------------------------------------------------

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
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock  (compiled into
    !     ssht library by JDM September 2005)
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
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock  (compiled into
    !     ssht library by JDM September 2005)
    !--------------------------------------------------------------------------

    subroutine ssht_dl_beta_operator_core(dl, beta, l)

      integer, intent(in) :: l
      real(kind = dp), intent(out) :: dl(-l:l,-l:l)
      real(kind = dp), intent(in) :: beta

  
      real(kind = dp) :: lambda, xllp1, xm, big, bigi, bigi2, c, s, omc, &
        lbig, expo, renorm, ratio, c2, t, si, ang, big_const
      real(kind = dp) :: cp(1:2 * l + 1)
      real(kind = dp) :: cpi(1:2 * l + 1)
      real(kind = dp) :: cp2(1:2 * l + 1)
      real(kind = dp) :: log_first_row(1:2 * l + 1)
      real(kind = dp) :: sign(1:2 * l + 1)
      real(kind = dp) :: lrenorm(1:2 * l + 1)

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
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock  (compiled into
    !     ssht library by JDM September 2005)
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


    !--------------------------------------------------------------------------
    ! Risbo's method
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! ssht_dl_beta_risbo_full
    !
    !! Calculates the *full* lth plane of a d-matrix using Risbo's
    !! recursion method.  For l>1, require the dl plane to be computed
    !! already with values for l-1.
    !!
    !! Variables:
    !!   - dl: Dl matrix values computed for lth plane.
    !!   - beta: Beta euler angle to compute dl matrix for.
    !!   - l: Plane of dl matrix to compute values for.
    !
    !! @author D. J. Mortlock
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock (compiled into
    !     ssht library by JDM December 2010)
    !--------------------------------------------------------------------------

    subroutine ssht_dl_beta_risbo_full(dl, beta, l)

      integer, intent(in) :: l
      real(kind = dp), intent(out) :: dl(-l:l,-l:l)
      real(kind = dp), intent(in) :: beta

      integer :: sign, i, j, k, m, mm
      real(kind = dp) :: sinb, cosb, sinhb, coshb, rj, dlj, ddj
      real(kind = dp) :: sqrt_jmk, sqrt_kp1, sqrt_jmi, sqrt_ip1
      real(kind=dp) :: dd(0: 2 * l + 1, 0: 2 * l + 1)

      if (l < 0) then

         call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_dl_beta_recursion', &
              comment_add='el < 0');

      ! For the special cases of l = 0 use the direct formula.
      else if (l == 0) then

         dl(0, 0) = 1.0_dp

      ! For the special cases of l = 1 use the direct formula.
      else if (l == 1) then

         ! These formulae are taken directly from Brink & Satchler.

         cosb = cos(beta)
         sinb = sin(beta)

         coshb = cos(beta / 2.0)
         sinhb = sin(beta / 2.0)

         dl(- 1, - 1) = coshb**2
         dl(- 1, 0) = sinb / SQRT2
         dl(- 1, 1) = sinhb**2

         dl(0, - 1) = - sinb / SQRT2
         dl(0, 0) = cosb
         dl(0, 1) = sinb / SQRT2

         dl(1, - 1) = sinhb**2
         dl(1, 0) = - sinb / SQRT2
         dl(1, 1) = coshb**2

      else

         sinhb = sin(beta / 2.0) ! p
         coshb = - cos(beta / 2.0) !-q

         ! Initialise the plane of the dl-matrix to 0.0 for the recursion
         ! from l - 1 to l - 1/2.

         dd(0: 2 * l + 1, 0: 2 * l + 1) = 0.0_dp

         j = 2 * l - 1
         rj = real(j)
         do k = 0, j - 1
            sqrt_jmk = dsqrt(real(j-k,kind=dp))
            sqrt_kp1 = dsqrt(real(k+1,kind=dp))
            do i = 0, j - 1
               sqrt_jmi = dsqrt(real(j-i,kind=dp))
               sqrt_ip1 = dsqrt(real(i+1,kind=dp))
               dlj = dl(k - (l - 1), i - (l - 1)) / rj
               dd(i, k) = dd(i, k) &
                    + sqrt_jmi * sqrt_jmk * dlj * coshb
               dd(i + 1, k) = dd(i + 1, k) &
                    - sqrt_ip1 * sqrt_jmk * dlj * sinhb
               dd(i, k + 1) = dd(i, k + 1) &
                    + sqrt_jmi * sqrt_kp1 * dlj * sinhb
               dd(i + 1, k + 1) = dd(i + 1, k + 1) &
                    + sqrt_ip1 * sqrt_kp1 * dlj * coshb
            end do
         end do

         ! Having constructed the d^(l+1/2) matrix in dd, do the second
         ! half-step recursion from dd to dl. Start by initilalising  
         ! the plane of the dl-matrix to 0.0.

         dl(- l: l, - l: l) = 0.0_dp

         j = 2 * l
         rj = real(j)
         do k = 0, j - 1
            sqrt_jmk = dsqrt(real(j-k,kind=dp))
            sqrt_kp1 = dsqrt(real(k+1,kind=dp))
            do i = 0, j - 1
               sqrt_jmi = dsqrt(real(j-i,kind=dp))
               sqrt_ip1 = dsqrt(real(i+1,kind=dp))
               ddj = dd(i, k) / rj
               dl(k - l, i - l) = dl(k - l, i - l) &
                    + sqrt_jmi * sqrt_jmk * ddj * coshb
               dl(k - l, i + 1 - l) = dl(k - l, i + 1 - l) &
                    - sqrt_ip1 * sqrt_jmk * ddj * sinhb
               dl(k + 1 - l, i - l) = dl(k + 1 - l, i - l) &
                    + sqrt_jmi * sqrt_kp1 * ddj * sinhb
               dl(k + 1 - l, i + 1 - l) = dl(k + 1 - l, i + 1 - l) &
                    + sqrt_ip1 * sqrt_kp1 * ddj * coshb
            end do
         end do

      end if

    end subroutine ssht_dl_beta_risbo_full


    !--------------------------------------------------------------------------
    ! ssht_dl_beta_risbo_full_table
    !
    !! Calculates the *full* lth plane of a d-matrix using Risbo's
    !! recursion method.  For l>1, require the dl plane to be computed
    !! already with values for l-1.  Also takes a table of precomputed
    !! square roots of integers to avoid recomputing them.
    !!
    !! Variables:
    !!   - dl: Dl matrix values computed for lth plane.
    !!   - beta: Beta euler angle to compute dl matrix for.
    !!   - l: Plane of dl matrix to compute values for.
    !!   - sqrt_tbl(0:2*l): Table of precomputed square roots of integers 
    !!     0:2*l. 
    !
    !! @author D. J. Mortlock
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock (compiled into
    !     ssht library by JDM December 2010)
    !--------------------------------------------------------------------------

    subroutine ssht_dl_beta_risbo_full_table(dl, beta, l, sqrt_tbl)

      integer, intent(in) :: l
      real(kind = dp), intent(out) :: dl(-l:l,-l:l)
      real(kind = dp), intent(in) :: beta
      real(kind = dp), intent(in) :: sqrt_tbl(0:2*l)

      integer :: sign, i, j, k, m, mm
      real(kind = dp) :: sinb, cosb, sinhb, coshb, rj, dlj, ddj
      real(kind = dp) :: sqrt_jmk, sqrt_kp1, sqrt_jmi, sqrt_ip1
      real(kind=dp) :: dd(0: 2 * l + 1, 0: 2 * l + 1)


      if (l < 0) then

         call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_dl_beta_recursion', &
              comment_add='el < 0');

      ! For the special cases of l = 0 use the direct formula.
      else if (l == 0) then

         dl(0, 0) = 1.0_dp

      ! For the special cases of l = 1 use the direct formula.
      else if (l == 1) then

         ! These formulae are taken directly from Brink & Satchler.

         cosb = cos(beta)
         sinb = sin(beta)

         coshb = cos(beta / 2.0)
         sinhb = sin(beta / 2.0)

         dl(- 1, - 1) = coshb**2
         dl(- 1, 0) = sinb / SQRT2
         dl(- 1, 1) = sinhb**2

         dl(0, - 1) = - sinb / SQRT2
         dl(0, 0) = cosb
         dl(0, 1) = sinb / SQRT2

         dl(1, - 1) = sinhb**2
         dl(1, 0) = - sinb / SQRT2
         dl(1, 1) = coshb**2

      else

         sinhb = sin(beta / 2.0) ! p
         coshb = - cos(beta / 2.0) !-q

         ! Initialise the plane of the dl-matrix to 0.0 for the recursion
         ! from l - 1 to l - 1/2.

         dd(0: 2 * l + 1, 0: 2 * l + 1) = 0.0_dp

         j = 2 * l - 1
         rj = real(j)
         do k = 0, j - 1
            do i = 0, j - 1
               dlj = dl(k - (l - 1), i - (l - 1)) / rj
               dd(i, k) = dd(i, k) &
                    + sqrt_tbl(j-i) * sqrt_tbl(j-k) * dlj * coshb
               dd(i + 1, k) = dd(i + 1, k) &
                    - sqrt_tbl(i+1) * sqrt_tbl(j-k) * dlj * sinhb
               dd(i, k + 1) = dd(i, k + 1) &
                    + sqrt_tbl(j-i) * sqrt_tbl(k+1) * dlj * sinhb
               dd(i + 1, k + 1) = dd(i + 1, k + 1) &
                    + sqrt_tbl(i+1) * sqrt_tbl(k+1) * dlj * coshb
            end do
         end do

         ! Having constructed the d^(l+1/2) matrix in dd, do the second
         ! half-step recursion from dd to dl. Start by initilalising  
         ! the plane of the dl-matrix to 0.0.

         dl(- l: l, - l: l) = 0.0_dp

         j = 2 * l
         rj = real(j)
         do k = 0, j - 1
            do i = 0, j - 1
               ddj = dd(i, k) / rj
               dl(k - l, i - l) = dl(k - l, i - l) &
                    + sqrt_tbl(j-i) * sqrt_tbl(j-k) * ddj * coshb
               dl(k - l, i + 1 - l) = dl(k - l, i + 1 - l) &
                    - sqrt_tbl(i+1) * sqrt_tbl(j-k) * ddj * sinhb
               dl(k + 1 - l, i - l) = dl(k + 1 - l, i - l) &
                    + sqrt_tbl(j-i) * sqrt_tbl(k+1) * ddj * sinhb
               dl(k + 1 - l, i + 1 - l) = dl(k + 1 - l, i + 1 - l) &
                    + sqrt_tbl(i+1) * sqrt_tbl(k+1) * ddj * coshb
            end do
         end do

      end if

    end subroutine ssht_dl_beta_risbo_full_table


    !--------------------------------------------------------------------------
    ! ssht_dl_beta_risbo_half_table
    !
    !! Calculates *half* (for m = -l:l and mm = -l:0) of lth plane of
    !! a d-matrix using Risbo's recursion method.  For l>1, require
    !! the dl plane to be computed already with values for l-1.  Also
    !! takes a table of precomputed square roots of integers to avoid
    !! recomputing them.
    !!
    !! Variables:
    !!   - dl: Dl matrix values computed for lth plane.
    !!   - beta: Beta euler angle to compute dl matrix for.
    !!   - l: Plane of dl matrix to compute values for.
    !!   - sqrt_tbl(0:2*l): Table of precomputed square roots of integers 
    !!     0:2*l. 
    !
    !! @author D. J. Mortlock
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock (compiled into
    !     ssht library by JDM December 2010)
    !--------------------------------------------------------------------------

    subroutine ssht_dl_beta_risbo_half_table(dl, beta, l, sqrt_tbl)

      integer, intent(in) :: l
      real(kind = dp), intent(out) :: dl(-l:l,-l:l)
      real(kind = dp), intent(in) :: beta
      real(kind = dp), intent(in) :: sqrt_tbl(0:2*l)

      integer :: sign, i, j, k, m, mm
      real(kind = dp) :: sinb, cosb, sinhb, coshb, rj, dlj, ddj
      real(kind = dp) :: sqrt_jmk, sqrt_kp1, sqrt_jmi, sqrt_ip1
      real(kind=dp) :: dd(0: 2 * l + 1, 0: 2 * l + 1)

      if (l < 0) then

         call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_dl_beta_recursion', &
              comment_add='el < 0');

      ! For the special cases of l = 0 use the direct formula.
      else if (l == 0) then

         dl(0, 0) = 1.0_dp

      ! For the special cases of l = 1 use the direct formula.
      else if (l == 1) then

         ! These formulae are taken directly from Brink & Satchler.

         cosb = cos(beta)
         sinb = sin(beta)

         coshb = cos(beta / 2.0)
         sinhb = sin(beta / 2.0)

         dl(- 1, - 1) = coshb**2
         dl(- 1, 0) = sinb / SQRT2
         dl(- 1, 1) = sinhb**2

         dl(0, - 1) = - sinb / SQRT2
         dl(0, 0) = cosb
         dl(0, 1) = sinb / SQRT2

         dl(1, - 1) = sinhb**2
         dl(1, 0) = - sinb / SQRT2
         dl(1, 1) = coshb**2

      else

         sinhb = sin(beta / 2.0) ! p
         coshb = - cos(beta / 2.0) !-q

         ! Initialise the plane of the dl-matrix to 0.0 for the recursion
         ! from l - 1 to l - 1/2.

         dd(0: 2 * l + 1, 0: 2 * l + 1) = 0.0_dp

         j = 2 * l - 1
         rj = real(j)
         do k = 0, j - 1
            do i = 0, l !j/2 + 1
               dlj = dl(k - (l - 1), i - (l - 1)) / rj
               dd(i, k) = dd(i, k) &
                    + sqrt_tbl(j-i) * sqrt_tbl(j-k) * dlj * coshb
               dd(i + 1, k) = dd(i + 1, k) &
                    - sqrt_tbl(i+1) * sqrt_tbl(j-k) * dlj * sinhb
               dd(i, k + 1) = dd(i, k + 1) &
                    + sqrt_tbl(j-i) * sqrt_tbl(k+1) * dlj * sinhb
               dd(i + 1, k + 1) = dd(i + 1, k + 1) &
                    + sqrt_tbl(i+1) * sqrt_tbl(k+1) * dlj * coshb
            end do
         end do

         ! Having constructed the d^(l+1/2) matrix in dd, do the second
         ! half-step recursion from dd to dl. Start by initilalising  
         ! the plane of the dl-matrix to 0.0.

         dl(- l: l, - l: l) = 0.0_dp

         j = 2 * l
         rj = real(j)
         do k = 0, j - 1
            do i = 0, l !j/2
               ddj = dd(i, k) / rj
               dl(k - l, i - l) = dl(k - l, i - l) &
                    + sqrt_tbl(j-i) * sqrt_tbl(j-k) * ddj * coshb
               dl(k - l, i + 1 - l) = dl(k - l, i + 1 - l) &
                    - sqrt_tbl(i+1) * sqrt_tbl(j-k) * ddj * sinhb
               dl(k + 1 - l, i - l) = dl(k + 1 - l, i - l) &
                    + sqrt_tbl(j-i) * sqrt_tbl(k+1) * ddj * sinhb
               dl(k + 1 - l, i + 1 - l) = dl(k + 1 - l, i + 1 - l) &
                    + sqrt_tbl(i+1) * sqrt_tbl(k+1) * ddj * coshb
            end do
         end do

      end if

    end subroutine ssht_dl_beta_risbo_half_table


    !--------------------------------------------------------------------------
    ! ssht_dl_beta_risbo_fill_half2full
    !
    !! Fill in the full Wigner plane from the half m = -l:l and 
    !! mm = -l:0, i.e. compute positive mm's.
    !!
    !! Variables:
    !!   - dl: Dl matrix values for lth plane.
    !!   - l: Plane of dl matrix to compute values for.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   December 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_dl_beta_risbo_fill_half2full(dl, l)

      integer, intent(in) :: l
      real(kind = dp), intent(inout) :: dl(-l:l,-l:l)

      integer :: m, mm

      ! Symmetry through origin.
      do m = -l, l
         do mm = 1, l
            dl(m,mm) = (-1)**(m+mm) * dl(-m,-mm)
         end do
      end do

    end subroutine ssht_dl_beta_risbo_fill_half2full


    !--------------------------------------------------------------------------
    ! ssht_dl_halfpi_risbo_eighth_table
    !
    !! Calculates *eighth* of lth plane of a d-matrix for PI/2 using Risbo's
    !! recursion method.  For l>1, require the dl plane to be computed
    !! already with values for l-1.  Also takes a table of precomputed
    !! square roots of integers to avoid recomputing them.
    !!
    !! Notes:
    !!   -  **** COMPUTES HALF PLANE AT PRESENT (NOT EIGHTH) ****
    !!
    !! Variables:
    !!   - dl: Dl matrix values computed for lth plane.
    !!   - beta: Beta euler angle to compute dl matrix for.
    !!   - l: Plane of dl matrix to compute values for.
    !!   - sqrt_tbl(0:2*l): Table of precomputed square roots of integers 
    !!     0:2*l. 
    !
    !! @author D. J. Mortlock
    !
    ! Revisions:
    !   September 2005 - Written by Daniel Mortlock (compiled into
    !     ssht library by JDM December 2010)
    !--------------------------------------------------------------------------

    subroutine ssht_dl_halfpi_risbo_eighth_table(dl, l, sqrt_tbl)

      integer, intent(in) :: l
      real(kind = dp), intent(out) :: dl(-l:l,-l:l)
      real(kind = dp), intent(in) :: sqrt_tbl(0:2*l)

      integer :: sign, i, j, k, m, mm
      real(kind = dp) :: sinb, cosb, sinhb, coshb, rj, dlj, ddj
      real(kind = dp) :: sqrt_jmk, sqrt_kp1, sqrt_jmi, sqrt_ip1
      real(kind=dp) :: dd(0: 2 * l + 1, 0: 2 * l + 1)
      real(kind = dp), parameter :: beta = PION2

      if (l < 0) then

         call ssht_error(SSHT_ERROR_ARG_INVALID, 'ssht_dl_beta_recursion', &
              comment_add='el < 0');

      ! For the special cases of l = 0 use the direct formula.
      else if (l == 0) then

         dl(0, 0) = 1.0_dp

      ! For the special cases of l = 1 use the direct formula.
      else if (l == 1) then

         ! These formulae are taken directly from Brink & Satchler.

         cosb = cos(beta)
         sinb = sin(beta)

         coshb = cos(beta / 2.0)
         sinhb = sin(beta / 2.0)

         dl(- 1, - 1) = coshb**2
         dl(- 1, 0) = sinb / SQRT2
         dl(- 1, 1) = sinhb**2

         dl(0, - 1) = - sinb / SQRT2
         dl(0, 0) = cosb
         dl(0, 1) = sinb / SQRT2

         dl(1, - 1) = sinhb**2
         dl(1, 0) = - sinb / SQRT2
         dl(1, 1) = coshb**2

      else

         sinhb = sin(beta / 2.0) ! p
         coshb = - cos(beta / 2.0) !-q

         ! Initialise the plane of the dl-matrix to 0.0 for the recursion
         ! from l - 1 to l - 1/2.

         dd(0: 2 * l + 1, 0: 2 * l + 1) = 0.0_dp

         j = 2 * l - 1
         rj = real(j)
!!$         do k = 0, j-1
         do k = 0, j - 1
!!$            do i = 0, j - 1
            do i = 0, l!j/2 + 1
               dlj = dl(k - (l - 1), i - (l - 1)) / rj
               dd(i, k) = dd(i, k) &
                    + sqrt_tbl(j-i) * sqrt_tbl(j-k) * dlj * coshb
               dd(i + 1, k) = dd(i + 1, k) &
                    - sqrt_tbl(i+1) * sqrt_tbl(j-k) * dlj * sinhb
               dd(i, k + 1) = dd(i, k + 1) &
                    + sqrt_tbl(j-i) * sqrt_tbl(k+1) * dlj * sinhb
               dd(i + 1, k + 1) = dd(i + 1, k + 1) &
                    + sqrt_tbl(i+1) * sqrt_tbl(k+1) * dlj * coshb
            end do
         end do

         ! Having constructed the d^(l+1/2) matrix in dd, do the second
         ! half-step recursion from dd to dl. Start by initilalising  
         ! the plane of the dl-matrix to 0.0.

         dl(- l: l, - l: l) = 0.0_dp

         j = 2 * l
         rj = real(j)
!!$         do k = 0, j-1
         do k = 0, j - 1
!!$            do i = 0, j - 1
            do i = 0, l!j/2
               ddj = dd(i, k) / rj
               dl(k - l, i - l) = dl(k - l, i - l) &
                    + sqrt_tbl(j-i) * sqrt_tbl(j-k) * ddj * coshb
               dl(k - l, i + 1 - l) = dl(k - l, i + 1 - l) &
                    - sqrt_tbl(i+1) * sqrt_tbl(j-k) * ddj * sinhb
               dl(k + 1 - l, i - l) = dl(k + 1 - l, i - l) &
                    + sqrt_tbl(j-i) * sqrt_tbl(k+1) * ddj * sinhb
               dl(k + 1 - l, i + 1 - l) = dl(k + 1 - l, i + 1 - l) &
                    + sqrt_tbl(i+1) * sqrt_tbl(k+1) * ddj * coshb
            end do
         end do

      end if

    end subroutine ssht_dl_halfpi_risbo_eighth_table


    !--------------------------------------------------------------------------
    ! ssht_dl_halfpi_risbo_fill_eighth2full
    !
    !! Fill in the full Wigner plane from eighth...
    !!
    !! Notes:
    !!   -  **** FILLS FROM HALF PLANE AT PRESENT (NOT EIGHTH) ****
    !!
    !! Variables:
    !!   - dl: Dl matrix values for lth plane.
    !!   - l: Plane of dl matrix to compute values for.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   December 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_dl_halfpi_risbo_fill_eighth2full(dl, l)

      integer, intent(in) :: l
      real(kind = dp), intent(inout) :: dl(-l:l,-l:l)

      integer :: m, mm

!!$      ! Symmetry in m.
!!$      do m = 1, l
!!$         do mm = -l, -1
!!$            dl(m,mm) = (-1)**(l+mm) * dl(-m,mm)
!!$         end do
!!$      end do
      ! Symmetry in mm.
      do m = -l, l
         do mm = 1, l
            dl(m,mm) = (-1)**(l+m) * dl(m,-mm)
         end do
      end do

    end subroutine ssht_dl_halfpi_risbo_fill_eighth2full


    !--------------------------------------------------------------------------
    ! Trapani & Navaza's method
    !--------------------------------------------------------------------------    
    
    !--------------------------------------------------------------------------
    ! ssht_dl_halfpi_trapani_eighth
    !
    !! Calculates *eighth* (for m = 0:l and mm = 0:m of lth plane of a
    !! d-matrix for PI/2 using Trapani & Navaza's recursion method.  For l>0,
    !! require the dl plane to be computed already with values for
    !! l-1.
    !!
    !! Variables:
    !!   - dl: Dl matrix values computed for lth plane.
    !!   - l: Plane of dl matrix to compute values for.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   December 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_dl_halfpi_trapani_eighth(dl, el)

      integer, intent(in) :: el
      real(dp), intent(inout) :: dl(-el:el,-el:el)

      real(dp) :: dmm(0:el)
      integer :: m, mm
      real(dp) :: t1, t2

      if (el == 0) then

         dl(0,0) = 1d0

      else

         ! Eqn (9) of T&N (2006).
         dmm(0) = - sqrt( (2d0*el-1d0) / real(2d0*el,dp) ) &
              * dl(el-1,0)

         ! Eqn (10) of T&N (2006).
         do mm = 1, el
            dmm(mm) = sqrt( el/2d0 * (2d0*el-1d0) / real((el+mm) * (el+mm-1), dp) ) &
                 * dl(el-1,mm-1)
         end do

         dl(el,0:el) = dmm(0:el)

         ! Eqn (11) of T&N (2006).
         do mm = 0, el

            ! m = el - 1 case (t2 = 0). 
            m = el-1
            t1 = ( 2e0 * mm / sqrt(real((el-m) * (el+m+1), dp)) ) * &
                 dl(m+1,mm)
            dl(m,mm) = t1

            ! Remaining m cases.
            do m = el-2, mm, -1
               t1 = ( 2e0 * mm / sqrt(real((el-m) * (el+m+1), dp)) ) * &
                    dl(m+1,mm)
               t2 = sqrt( (el-m-1) * (el+m+2) / real((el-m) * (el+m+1), dp) ) * &
                    dl(m+2,mm)
               dl(m,mm) = t1 - t2
            end do
         end do

      end if

    end subroutine ssht_dl_halfpi_trapani_eighth


    !--------------------------------------------------------------------------
    ! ssht_dl_halfpi_trapani_eighth_table
    !
    !! Calculates *eighth* (for m = 0:l and mm = 0:m of lth plane of a
    !! d-matrix for PI/2 using Trapani & Navaza's recursion method.  For l>0,
    !! require the dl plane to be computed already with values for
    !! l-1.  Also takes a table of precomputed square roots of
    !! integers to avoid recomputing them.
    !!
    !! Variables:
    !!   - dl: Dl matrix values computed for lth plane.
    !!   - l: Plane of dl matrix to compute values for.
    !!   - sqrt_tbl(0:2*l+1): Table of precomputed square roots of integers 
    !!     0:2*l+1. 
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   December 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_dl_halfpi_trapani_eighth_table(dl, el, sqrt_tbl)

      integer, intent(in) :: el
      real(dp), intent(inout) :: dl(0:el,0:el)
      real(kind = dp), intent(in) :: sqrt_tbl(0:2*el+1)


      real(dp) :: dmm(0:el)
      integer :: m, mm
      real(dp) :: t1, t2

      if (el == 0) then

         dl(0,0) = 1d0

      else

         ! Eqn (9) of T&N (2006).
         dmm(0) = - sqrt_tbl(2*el-1) / sqrt_tbl(2*el) &
              * dl(el-1,0)

         ! Eqn (10) of T&N (2006).
         do mm = 1, el
            dmm(mm) = sqrt_tbl(el) / SQRT2 &
                 * sqrt_tbl(2*el-1) / sqrt_tbl(el+mm) / sqrt_tbl(el+mm-1) &
                 * dl(el-1,mm-1)
         end do

         dl(el,0:el) = dmm(0:el)

         ! Eqn (11) of T&N (2006).
         do mm = 0, el

            ! m = el-1 case (t2 = 0). 
            m = el-1
            dl(m,mm) = 2e0 * mm / sqrt_tbl(el-m) / sqrt_tbl(el+m+1) &
                 * dl(m+1,mm)

            ! Remaining m cases.
            do m = el-2, mm, -1
               t1 = 2e0 * mm / sqrt_tbl(el-m) / sqrt_tbl(el+m+1) &
                    * dl(m+1,mm)
               t2 = sqrt_tbl(el-m-1) * sqrt_tbl(el+m+2) / sqrt_tbl(el-m) / sqrt_tbl(el+m+1) &
                    * dl(m+2,mm)
               dl(m,mm) = t1 - t2
            end do
         end do

      end if

    end subroutine ssht_dl_halfpi_trapani_eighth_table


    !--------------------------------------------------------------------------
    ! ssht_dl_halfpi_trapani_fill_eighth2full
    !
    !! Fill in the full Wigner plane from the eighth m = 0:l and 
    !! mm = 0:m.
    !!
    !! Variables:
    !!   - dl: Dl matrix values for lth plane.
    !!   - l: Plane of dl matrix to compute values for.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   December 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------
 
    subroutine ssht_dl_halfpi_trapani_fill_eighth2full(dl, el)

      integer, intent(in) :: el
      real(kind = dp), intent(inout) :: dl(-el:el,-el:el)

      integer :: m, mm

      ! Diagonal symmetry to fill in quarter.
      do m = 0, el
         do mm = m+1, el
            dl(m,mm) = (-1)**(m+mm) * dl(mm,m)
         end do
      end do

      ! Symmetry in m to fill in half.
      do m = -el, -1
         do mm = 0, el
            dl(m,mm) = (-1)**(el+mm) * dl(-m,mm)
         end do
      end do

      ! Symmetry in mm to fill in remaining plane.
      do m = -el, el
         do mm = -el, -1
            dl(m,mm) = (-1)**(el+m) * dl(m,-mm)
         end do
      end do

    end subroutine ssht_dl_halfpi_trapani_fill_eighth2full


    !--------------------------------------------------------------------------
    ! ssht_dl_halfpi_trapani_fill_eighth2quarter
    !
    !! Fill in quarter Wigner plane for m = 0:l and mm = 0:l from the
    !! eighth m = 0:l and mm = 0:m, i.e. compute positive mm's.
    !!
    !! Variables:
    !!   - dl: Dl matrix values for lth plane.
    !!   - l: Plane of dl matrix to compute values for.
    !
    !! @author J. D. McEwen
    !
    ! Revisions:
    !   December 2010 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_dl_halfpi_trapani_fill_eighth2quarter(dl, el)

      integer, intent(in) :: el
!!$      real(kind = dp), intent(inout) :: dl(-el:el,-el:el)
      real(kind = dp), intent(inout) :: dl(0:el,0:el)

      integer :: m, mm
      real(dp) :: signs(0:el+1)

      ! Compute signs.
      do m = 0, el, 2
         signs(m)   =  1.0_dp
         signs(m+1) = -1.0_dp
      enddo

      ! Diagonal symmetry to fill in quarter.
      do m = 0, el
         do mm = m+1, el
            dl(m,mm) = signs(m) * signs(mm) * dl(mm,m)
         end do
      end do
!!$      do mm = 0, el
!!$         do m = 0, mm-1
!!$            dl(m,mm) = signs(m) * signs(mm) * dl(mm,m)
!!$         end do
!!$      end do

    end subroutine ssht_dl_halfpi_trapani_fill_eighth2quarter



!eighth2half



end module ssht_dl_mod
