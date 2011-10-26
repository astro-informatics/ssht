! SSHT package to perform spin spherical harmonic transforms
! Copyright (C) 2011  Jason McEwen
! See LICENSE.txt for license details


!------------------------------------------------------------------------------
! ssht_about
!
!> Prints information about the SSHT package, including version
!! and build numbers. 
!!
!! Usage: ssht_about
!!     
!! @author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!
! Revisions:
!   October 2011 - Written by Jason McEwen 
!------------------------------------------------------------------------------

program ssht_about

  implicit none

  write(*,'(a)') "=========================================================="
  write(*,'(a)') "SSHT package to perform spin spherical harmonic transforms"
  write(*,'(a)') "By Jason McEwen and Yves Wiaux"

  write(*,'(a)') "See www.jasonmcewen.org for more information."
  write(*,'(a)') "See LICENSE.txt for license details."

  write(*,'(a,a)') "Version: ", SSHT_VERSION
  write(*,'(a,a)') "Build: ", SSHT_BUILD
  write(*,'(a)') "=========================================================="

end program ssht_about

