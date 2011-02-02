! D-matrix for rotating spherical harmonics.

module general_dl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use general_types
  use general_const
  use general_error
  use general_maths

  implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Minimum Euler angle for which a rotation is calculated.

  real(kind = gndp), parameter :: GNEPS_EULER = 1.0e-5_gndp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculate the entire d^l_mm matrix.

  function gn_d_beta(beta, lmax) result(d)
    real(kind = gndp), intent(in) :: beta
    integer, intent(in) :: lmax
    real(kind = gndp), pointer :: d(:, :, :), dl(:, :)

    integer :: l, mm, m

    ! Allocate and initialise d-matrix and dl-matrix.

    allocate(dl(- lmax: lmax, - lmax: lmax))
    dl(- lmax: lmax, - lmax: lmax) = 0.0_gndp

    allocate(d(0: lmax, - lmax: lmax, - lmax: lmax))
    d(0: lmax, - lmax: lmax, - lmax: lmax) = 0.0_gndp

    ! Loop over the planes of the matrix, calculating each in turn.

    do l = 0, lmax

      call gn_dl_beta(dl, beta, l, 'recursion')

      do m = - l, l
        do mm = - l, l
          d(l, m, mm) = dl(m, mm)
        end do
      end do

    end do

    deallocate(dl)
    
  end function gn_d_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Tests the orthogonality of the columns or rows of the lth plane of 
  ! a d-matrix. The return value is the magnitude of the greatest of 
  ! ``dot-product'' of any two rows or columns.

  function gn_dl_orthtest(dl, l, rowcol) result(dotmax)
    real(kind = gndp), pointer :: dl(:, :)
    integer, intent(in) :: l
    character(len = *), intent(in) :: rowcol
    real(kind = gndp) :: dotmax

    real(kind = gndp) :: dot
    integer :: m, m1, m2, mm, mm1, mm2

    dotmax = 0.0_gndp

    if (rowcol == 'row') then

      do m1 = - l, l
        do m2 = m1 + 1, l  

          ! Work out the dot product of row_1 and row_2.

          dot = 0.0_gndp
          do mm = - l, l
            dot = dot + dl(m1, mm) * dl(m2, mm)
          end do   
          
          if (abs(dot) > dotmax) then
            dotmax = abs(dot)
          end if

        end do
      end do

    else if (rowcol == 'col') then

      do mm1 = - l, l
        do mm2 = mm1 + 1, l  

          ! Work out the dot product of col_1 and col_2.

          dot = 0.0_gndp   
          do m = - l, l
            dot = dot + dl(m, mm1) * dl(m, mm2)
          end do   
          
          if (abs(dot) > dotmax) then
            dotmax = abs(dot)
          end if

        end do
      end do

    else
      call gn_fatal('gn_dl_orthtest: rowcol unknown', rowcol)
    end if

  end function gn_dl_orthtest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Tests the normality of the columns or rows of the lth plane of a 
  ! d-matrix. The return value is the magnitude of the greatest of 
  ! ``dot-product'' of any row or column with itself.

  function gn_dl_normtest(dl, l, rowcol) result(dotmax)
   real(kind = gndp), pointer :: dl(:, :)
    integer, intent(in) :: l
    character(len = *), intent(in) :: rowcol
    real(kind = gndp) :: dotmax

    real(kind = gndp) :: dot
    integer :: m, mm

    dotmax = 0.0

    if (rowcol == 'row') then

      do m = - l, l

        ! Work out the dot product of row and row.

        dot = 0.0_gndp
        do mm = - l, l
          dot = dot + dl(m, mm)**2
        end do
         
        if (abs(dot - 1.0_gndp) > dotmax) then
          dotmax = abs(dot - 1.0_gndp)
        end if

      end do

    else if (rowcol == 'col') then

      do mm = - l, l

        ! Work out the dot product of col_1 and col_2.

        dot = 0.0_gndp  
        do m = - l, l
          dot = dot + dl(m, mm)**2
        end do
         
        if (abs(dot - 1.0_gndp) > dotmax) then
          dotmax = abs(dot - 1.0_gndp)
        end if

      end do

    else
      call gn_fatal('gn_dl_normtest: rowcol unknown', rowcol)
    end if

  end function gn_dl_normtest

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Use gn_dl_beta recursively to generate the lth plane of the 
  ! d^l_mm matrix. The memory for dl is allocated here, and this
  ! is only useful for low-l applications where the planes can't
  ! be used successively, except if the direct evaluation method 
  ! is used, in which case this is viable.

  function gn_dl_beta_f(beta, l, method) result(dl)
    real(kind = gndp), intent(in) :: beta
    integer, intent(in) :: l
    real(kind = gndp), pointer :: dl(:, :)
    character(len = *), intent(in) :: method

    integer :: ll

    ! Allocate enough memory for the final d-matrix.

    allocate(dl(- l: l, - l: l))
    dl = 0.0_gndp

    ! Evaluate the matrices, depending on the method of calculation.

    if (method == 'recursion') then
      do ll = 0, l
        call gn_dl_beta(dl, beta, ll, 'recursion')
      end do
    else
      call gn_dl_beta(dl, beta, l, method)
    end if

  end function gn_dl_beta_f

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculate the lth plane of the d-matrix appropriate to rotating
  ! real harmonic coefficients, as opposed to the normal complex 
  ! coefficients. This, for the moment, is a wrapper to the complex
  ! d-matrices, calculating them first and then converting. Hence 
  ! it is not clear how best to handle the recursion case; for the
  ! moment we'll assume that dl comes in with the values for the 
  ! real (l - 1)th d-matrix and that the conversion will be done.
  ! At any rate the d-matrix is assumed to be allocated already,
  ! with both indices running from - l to l.

  subroutine gn_redl_beta(redl, beta, l, method)
    real(kind = gndp), pointer :: redl(:, :)
    real(kind = gndp), intent(in) :: beta
    integer, intent(in) :: l
    character(len = *), intent(in) :: method

    integer :: m, mm
    real(kind = gndp) :: mmsign
    real(kind = gndp), pointer :: dl(:, :) => null()

    if (abs(beta) < GNEPS_EULER) then

      ! In this case the d-matrix is just the identity.

      call gn_warning('gn_redl_beta: beta = 0.0')

      do m = - l, l
        do mm = - l, l
          if (m == mm) then
            redl(m, mm) = 1.0
          else 
            redl(m, mm) = 0.0
          end if
        end do
      end do

    else 

      ! First up calculate the complex d-matrix for the relevant l.

      if (method == 'recursion') then
        call gn_fatal('gn_redl_beta: recursion method not yet programmed')
      else
        dl => gn_dl_beta_f(beta, l, method)
      end if
  
      ! Then convert over to the real d-matrix, section by section, 
      ! first initialising it to zeros.
  
      do m = - l, l
        do mm = - l, l
          redl(m, mm) = 0.0_gndp
        end do
      end do
  
      do m = - 1, - l, - 1
        mmsign = 1.0
        do mm = - 1, - l, - 1
          mmsign = - mmsign
          redl(m, mm) = dl(- m, - mm) - mmsign * dl(- m, mm)
        end do
      end do
  
      m = 0
      mm = 0
      redl(0, 0) = dl(0, 0)
  
      m = 0
      mmsign = 1.0
      do mm = 1, l
        mmsign = - mmsign
        redl(0, mm) = (dl(0, mm) + mmsign * dl(0, - mm)) / GNSQRT2
      end do
  
      mm = 0
      do m = 1, l
        redl(m, 0) = GNSQRT2 * dl(m, 0)
      end do
  
      do m = 1, l
        mmsign = 1.0
        do mm = 1, l
          mmsign = - mmsign
          redl(m, mm) = dl(m, mm) + mmsign * dl(m, - mm)   
        end do
      end do
  
      ! Finally free up the temporary memory. 

      deallocate(dl)
 
    end if

  end subroutine gn_redl_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculate the lth plane of the d-matrix to rotate spherical harmonic
  ! coefficients about the (positive) y-axis. This is just a wrapper to 
  ! the specific methods of doing this: Risbo's recursion; Turok's 
  ! operator-based method; and the direct method (from Brink & Satchler).
  ! This assumes that dl is already allocated to have both indices 
  ! running at least from - l to l. Further, in the case of the recursion
  ! method, this requires that, on entry, dl contain the (l - 1)th 
  ! plane of the the d-matrix.

  subroutine gn_dl_beta(dl, beta, l, method)
    real(kind = gndp), pointer :: dl(:, :)
    real(kind = gndp), intent(in) :: beta
    integer, intent(in) :: l
    character(len = *), intent(in) :: method

    if (method == 'recursion') then

      ! Use Risbo's recursion (which requires that dl have the (l - 1)th
      ! d-matrix in place.

      call gn_dl_beta_recursion(dl, beta, l)

    else if (method == 'operator') then

      ! Use Turok & Bucher's operator-based method, which does not 
      ! require values on input. 

      call gn_dl_beta_operator(dl, beta, l)

    else if (method == 'direct') then

      ! Use a direct evaluation from Brink & Satchler.

      call gn_dl_beta_direct(dl, beta, l)

    else

      call gn_fatal('gn_dl_beta: method unknown', method)

    end if

  end subroutine gn_dl_beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculate the lth plane of a d-matrix using Turok & Bucher's 
  ! operator-based method. Note that the recursion used becomes 
  ! unstable at angles close to any integer multiple of pi (i.e.,
  ! when d^l_mm becomes close to diagonal and the rotation angle
  ! is minimal. Just when this happens hasn't been tested as yet, 
  ! although for beta = 10 deg and l = 400 everything seems fine,
  ! so a possibility is beta = 1 deg at l = 1000. When the answer
  ! to this question is known with some confidence, there will need
  ! to be some if statement to evaluate the matrix some other way,
  ! or at least warn the user that they're entering rocky ground.

  subroutine gn_dl_beta_operator(dl, beta, l)
    real(kind = gndp), pointer :: dl(:, :)
    real(kind = gndp), intent(in) :: beta
    integer, intent(in) :: l

    ! Fill in the top quarter of the d-matrix.

    call gn_dl_beta_operator_core(dl, beta, l)

    ! Use its symmetry properties to fill in the rest of it.

    call gn_dl_beta_operator_fill(dl, l)

  end subroutine gn_dl_beta_operator

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Does the left quarter of the d-matrix. Beta is the angle of rotation,
  ! l is the plane of the matrix required and dl is the two dimensional 
  ! array representing the lth plane of the matrix. This array must be
  ! already allocated on entry with dimensions given by the command
  ! allocate(dl(- lmax: lmax, - lmax: lmax)) where lmax >= l.

  subroutine gn_dl_beta_operator_core(dl, beta, l)
    real(kind = gndp), pointer :: dl(:, :)
    real(kind = gndp), intent(in) :: beta
    integer, intent(in) :: l

    real(kind = gndp) :: lambda, xllp1, xm, big, bigi, bigi2, c, s, omc, &
      lbig, expo, renorm, ratio, c2, t, si, ang, big_const
    real(kind = gndp), pointer :: cp(:) => null(), cpi(:) => null(), &
      cp2(:) => null(), log_first_row(:) => null(), sign(:) => null(), &
      lrenorm(:) => null()
    integer :: index, i, m, im, j, lp1

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

    log_first_row(1)=(2.d0*real(l,kind=gndp))*dlog(dabs(c2))
    sign(1)=1.d0
    do i=2, 2*l+1
      m=l+1-i
      ratio=dsqrt(real(l+m+1,kind=gndp)/real(l-m,kind=gndp))
      log_first_row(i)=log_first_row(i-1) &
        +dlog(ratio)+dlog(dabs(t))
      sign(i)=sign(i-1)*t/dabs(t)
    end do

    big=big_const
    lbig=dlog(big)
    bigi=1.d0/big_const
    bigi2=1.d0/big_const**2
    xllp1=real(l*(l+1),kind=gndp)

    ! Initialising coefficients cp(m)= cplus(l-m).

    do m=1,l
      xm=real(l-m,kind=gndp)
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
      lambda=(real(l+1,kind=gndp)*omc-real(index,kind=gndp)+c)*si
      dl(index - lp1, 2 - lp1)=lambda*dl(index - lp1, 1 - lp1)*cpi(1)
      if(index.gt.2) then
        do m=2,index-1
          lambda=(real(l+1,kind=gndp)*omc &
            -real(index,kind=gndp)+real(m,kind=gndp)*c)/s
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
      lambda=(real(l+1,kind=gndp)*omc-real(index,kind=gndp)+c)/s
      dl(index - lp1, 2 - lp1)=lambda*dl(index - lp1, 1 - lp1)*cpi(1)
      if(index.lt.2*l) then
        do m=2,2*l-index+1
          lambda=(real(l+1,kind=gndp)*omc-real(index,kind=gndp) &
            +real(m,kind=gndp)*c)/s
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

!!!!!
! Free memory.
! (Added by Jason McEwen 07/06/2004)
    deallocate(cp)
    deallocate(cpi)
    deallocate(cp2)
    deallocate(log_first_row)
    deallocate(sign)
    deallocate(lrenorm)
!!!!!!


  end subroutine gn_dl_beta_operator_core

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Computes the three remaining (top, bottom and right) quarters of the 
  ! d-matrix, given the left quarter.

  subroutine gn_dl_beta_operator_fill(dl, l)
    real(kind = gndp), pointer :: dl(:, :)
    integer, intent(in) :: l

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

  end subroutine gn_dl_beta_operator_fill

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculate the lth plane of a d-matrix using Brink & Satchler's 
  ! direct evaluation. This is only useful for testing purposes of 
  ! more efficient methods, and in fact breaks down at ls of between
  ! 10 and 100, depending on beta.

  subroutine gn_dl_beta_direct(dl, beta, l)
    real(kind = gndp), pointer :: dl(:, :)
    real(kind = gndp), intent(in) :: beta
    integer, intent(in) :: l

    real(kind = gndp) :: cosb, sinb, coshb, sinhb, tsign, rt, rl, rm, rmm, &
      logsinhb, logcoshb, csign, ssign
    integer :: t_min, t_max, t, m, mm, coshbsign, sinhbsign

    ! Test precomputed square roots and factorials.

!   call gn_sqt_test('gn_dl_beta_direct', 8)
!   call gn_logfact_test('gn_dl_beta_direct', 2 * l)

    if (l < 0) then

      ! No such thing as l < 0 d-matrix.

      call gn_fatal('gn_dl_beta_direct: l < 0', l)

    else if (l == 0) then

      ! The normalisation of a scalar function doesn't change upon rotation,
      ! hence d^0_m,m' is the identity.

      dl(0, 0) = 1.0

    else if (l == 1) then

      ! These formulae are taken directly from Brink & Satchler.

      cosb = cos(beta)
      sinb = sin(beta)

      coshb = cos(beta / 2.0)
      sinhb = sin(beta / 2.0)

      dl(- 1, - 1) = coshb**2
      dl(- 1, 0) = sinb / GNSQRT2
      dl(- 1, 1) = sinhb**2

      dl(0, - 1) = - sinb / GNSQRT2
      dl(0, 0) = cosb
      dl(0, 1) = sinb / GNSQRT2

      dl(1, - 1) = sinhb**2
      dl(1, 0) = - sinb / GNSQRT2
      dl(1, 1) = coshb**2

    else if (l == 2) then

      ! These formulae are taken directly from Brink & Satchler.

      cosb = cos(beta)
      sinb = sin(beta)

      coshb = cos(beta / 2.0)
      sinhb = sin(beta / 2.0)

      dl(- 2, - 2) = coshb**4
      dl(- 2, - 1) = 0.5 * sinb * (1 + cosb)
      dl(- 2, 0) = gn_sqt(3) / gn_sqt(8) * sinb**2
      dl(- 2, 1) = - 0.5 * sinb * (cosb - 1.0)
      dl(- 2, 2) = sinhb**4

      dl(- 1, - 2) = - 0.5 * sinb * (1 + cosb)
      dl(- 1, - 1) = 0.5 * (2.0 * cosb - 1.0) * (cosb + 1.0)
      dl(- 1, 0) = gn_sqt(3) / gn_sqt(2) * sinb * cosb
      dl(- 1, 1) = 0.5 * (2.0 * cosb + 1.0) * (1.0 - cosb)
      dl(- 1, 2) = - 0.5 * sinb * (cosb - 1.0)

      dl(0, - 2) = gn_sqt(3) / gn_sqt(8) * sinb**2
      dl(0, - 1) = - gn_sqt(3) / gn_sqt(2) * sinb * cosb
      dl(0, 0) = 0.5 * (3.0 * cosb**2 - 1.0)
      dl(0, 1) = gn_sqt(3) / gn_sqt(2) * sinb * cosb
      dl(0, 2) = gn_sqt(3) / gn_sqt(8) * sinb**2

      dl(1, - 2) = 0.5 * sinb * (cosb - 1.0)
      dl(1, - 1) = 0.5 * (2.0 * cosb + 1.0) * (1.0 - cosb)
      dl(1, 0) = - gn_sqt(3) / gn_sqt(2) * sinb * cosb
      dl(1, 1) = 0.5 * (2.0 * cosb - 1.0) * (cosb + 1.0)
      dl(1, 2) = 0.5 * sinb * (1 + cosb)

      dl(2, - 2) = sinhb**4
      dl(2, - 1) = 0.5 * sinb * (cosb - 1.0)
      dl(2, 0) = gn_sqt(3) / gn_sqt(8) * sinb**2
      dl(2, 1) = - 0.5 * sinb * (1 + cosb)
      dl(2, 2) = coshb**4

    else 

      ! Evaluate the direct summation provided by Brink & Satchler by
      ! using precomputed factorials. In this method each element is
      ! computed separately, and so the two outer loops simply cover
      ! the entire matrix.

      coshb = cos(beta / 2.0)
      if (coshb < 0.0) then
        coshbsign = - 1
      else 
        coshbsign = 1
      end if
      coshb = abs(coshb)
      logcoshb = log(abs(coshb))

      sinhb = sin(beta / 2.0)
      if (sinhb < 0.0) then
        sinhbsign = - 1
      else 
        sinhbsign = 1
      end if
      logsinhb = log(abs(sinhb))

      rl = real(l)

      do m = - l, l
        rm = real(m)

        do mm = - l, l
          rmm = real(mm)

          ! For each element set to zero and then evaluate Brink & 
          ! Satchler's sum over t, making sure that the factorials 
          ! are never negative.

          dl(m, mm) = 0.0

          t_min = max(0, m - mm)
          t_max = min(l + m, l - mm)

          if (gn_iseven(t_min)) then
            tsign = - 1.0
          else
            tsign = 1.0
          end if 

          do t = t_min, t_max
            rt = real(t)
            tsign = - tsign

            if ((coshbsign == - 1) &
              .and. (gn_isodd(2 * l + m - mm - 2 * t))) then
              csign = - 1.0
            else
              csign = 1.0
            end if

            if ((sinhbsign == - 1) &
              .and. (gn_isodd(2 * t + mm - m))) then
              ssign = - 1.0
            else
              ssign = 1.0
            end if

            dl(m, mm) = dl(m, mm) &
              + csign * ssign * tsign * exp(0.5 &
              * (gn_logfact(l + m) + gn_logfact(l - m) &
              + gn_logfact(l + mm) + gn_logfact(l - mm)) &
              - (gn_logfact(l + m - t) + gn_logfact(l - mm - t) &
              + gn_logfact(t) + gn_logfact(t + mm - m)) &
              + (2.0 * rl + rm - rmm - 2.0 * rt) * logcoshb &
              + (2.0 * rt + rmm - rm) * logsinhb) 

          end do

        end do

      end do 

    end if

  end subroutine gn_dl_beta_direct

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Calculate the lth plane of a rotational d^l_mm matrix using Risbo's
  ! recursion method. The result, which is d^l(m, m') in that order,
  ! is written onto the already-allocated dl_mm array, which, on entry,
  ! must contain dl-1_mm, from which the lth plane is calculated using
  ! a four-way  recursion ``through'' the l-1/2th plane (stored temporarily
  ! in array dd. Here beta is the polar angle of rotation (i.e., the 
  ! second Euler angle). This routine has been tested qualitatively, and 
  ! appears to rotate a field in the correct manner (i.e., when the field is 
  ! represented in real space, it is unchanged except in overall orientation)
  ! and the rows/columns of the the dl appear to be orthogonal and 
  ! normalised to unity. (Specifically, the error in orthogonality is 
  ! approximately delta = 10^(-17) * l for l up to several thousand, and the 
  ! relative error in the normalisation is slightly larger, being 
  ! roughly delta = 10^(-16) * l, also for l up to several thousand.
  ! Note also that dl must be allocated from - l to l in both dimensions
  ! on entry. Also, the order of the inner loop is unimportant numerically, 
  ! but, in FORTRAN, important for speed. The ordering as is is superior as 
  ! the inner loop is over the closest index.

  subroutine gn_dl_beta_recursion(dl, beta, l)
    real(kind = gndp), pointer :: dl(:, :)
    real(kind = gndp), intent(in) :: beta
    integer, intent(in) :: l

    integer :: sign, i, j, k, m, mm
    real(kind = gndp) :: sinhb, coshb, rj, dlj, ddj
    real(kind = gndp), pointer :: dd(:, :) => null()

    if (l < 0) then

      call gn_fatal('gn_dl_beta_recursion: l < 0', l)

    else if (l <= 1) then

      ! For the special cases of l = 0 and l = 1, use the direct formula.

      call gn_dl_beta_direct(dl, beta, l)

    else

      ! Check that the square root table has been initialised to sufficiently
      ! large values, and initialise the dd-matrix, used to store the
      ! half integer step. And precompute trigonometric functions.

!     call gn_sqt_test('gn_dl_beta_recursion', 2 * l + 1)
      allocate(dd(0: 2 * l + 1, 0: 2 * l + 1))
      sinhb = sin(beta / 2.0)
      coshb = - cos(beta / 2.0)

      ! Initialise the plane of the dl-matrix to 0.0 for the recursion
      ! from l - 1 to l - 1/2.

      dd(0: 2 * l + 1, 0: 2 * l + 1) = 0.0_gndp
 
      j = 2 * l - 1
      rj = real(j)
      do k = 0, j - 1
        do i = 0, j - 1
          dlj = dl(k - (l - 1), i - (l - 1)) / rj
          dd(i, k) = dd(i, k) &
            + gn_sqt(j - i) * gn_sqt(j - k) * dlj * coshb
          dd(i + 1, k) = dd(i + 1, k) &
            - gn_sqt(i + 1) * gn_sqt(j - k) * dlj * sinhb
          dd(i, k + 1) = dd(i, k + 1) &
            + gn_sqt(j - i) * gn_sqt(k + 1) * dlj * sinhb
          dd(i + 1, k + 1) = dd(i + 1, k + 1) &
            + gn_sqt(i + 1) * gn_sqt(k + 1) * dlj * coshb
        end do
      end do

      ! Having constructed the d^(l+1/2) matrix in dd, do the second
      ! half-step recursion from dd to dl. Start by initilalising  
      ! the plane of the dl-matrix to 0.0.

      dl(- l: l, - l: l) = 0.0_gndp
 
      j = 2 * l
      rj = real(j)
      do k = 0, j - 1
        do i = 0, j - 1
          ddj = dd(i, k) / rj
          dl(k - l, i - l) = dl(k - l, i - l) &
            + gn_sqt(j - i) * gn_sqt(j - k) * ddj * coshb
          dl(k - l, i + 1 - l) = dl(k - l, i + 1 - l) &
            - gn_sqt(i + 1) * gn_sqt(j - k) * ddj * sinhb
          dl(k + 1 - l, i - l) = dl(k + 1 - l, i - l) &
            + gn_sqt(j - i) * gn_sqt(k + 1) * ddj * sinhb
          dl(k + 1 - l, i + 1 - l) = dl(k + 1 - l, i + 1 - l) &
            + gn_sqt(i + 1) * gn_sqt(k + 1) * ddj * coshb
        end do
      end do
 
      ! Free up the temporary matrix.

      deallocate(dd)

    end if

  end subroutine gn_dl_beta_recursion

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module general_dl
