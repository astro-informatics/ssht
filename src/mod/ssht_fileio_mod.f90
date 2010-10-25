!------------------------------------------------------------------------------
! ssht_fileio_mod  -- SSHT library file IO class
! 
!! Functionality to read and write SSHT formatted fits and matlab files 
!! containing wavelet and scaling coefficients.  Both statically and 
!! dynamically allocated wavelet coefficients may be written and read from 
!! files (both data types have the same fits file format).
!
!! @author J. D. McEwen (mcewen@mrao.cam.ac.uk)
!! @version 0.1 November 2007
!
! Revisions:
!   November 2007 - Written by Jason McEwen
!------------------------------------------------------------------------------

module ssht_fileio_mod

  use ssht_types_mod
  use ssht_error_mod
  use ssht_core_mod

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: &
    ssht_fileio_matlab_wav_write, &
    ssht_fileio_matlab_wav_read, &
		ssht_fileio_fits_wav_write, &
		ssht_fileio_fits_wav_read


  !---------------------------------------
  ! Interfaces
  !---------------------------------------

  interface ssht_fileio_matlab_wav_write
     module procedure &
       ssht_fileio_matlab_wav_write_static, &
       ssht_fileio_matlab_wav_write_dynamic
  end interface

  interface ssht_fileio_matlab_wav_read
     module procedure &
       ssht_fileio_matlab_wav_read_static, &
       ssht_fileio_matlab_wav_read_dynamic
  end interface

  interface ssht_fileio_fits_wav_write
     module procedure &
       ssht_fileio_fits_wav_write_static, &
       ssht_fileio_fits_wav_write_dynamic
  end interface

  interface ssht_fileio_fits_wav_read
     module procedure &
       ssht_fileio_fits_wav_read_static, &
       ssht_fileio_fits_wav_read_dynamic
  end interface


  !----------------------------------------------------------------------------

  contains


    !--------------------------------------------------------------------------
    ! Matlab file IO 
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! ssht_fileio_matlab_wav_write_static
    !
    !! Writes (statically allocated) wavelet and scaling coefficients to an 
    !! output SSHT formatted .m matlab file and corresponding .dat data files. 
    !!
    !! Variables:
    !!  - wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1): Wavelet coefficients [input].
		!!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
		!!    [input].
		!!  - J: Maximum analysis scale depth [input].
		!!  - B: Harmonic band limit [input].
		!!  - N: Azimuthal band limit [input].
		!!  - bl_scoeff: Upper band limit for scaling coefficients [input].
		!!  - alpha: Basis dilation factor [input].
    !!  - filename: Name of the output matlab file to write [input].
    !!  - [comment]: Optional comment string to be added to the output matlab 
    !!    file header if present [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - May 2008
    !
    ! Revisions:
    !   May 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_matlab_wav_write_static(wav, scoeff, J, B, N, &
      bl_scoeff, alpha, filename, comment)

			integer, intent(in) :: J
			integer, intent(in) :: B
			integer, intent(in) :: N
			integer, intent(in) :: bl_scoeff
			real(dp), intent(in) :: alpha
			real(dp), intent(in) :: wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1)
			complex(dpc), intent(in) :: scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      integer :: fileid, fileid_wav, fileid_scoeff
      integer :: el, m, jj, aa, bb, gg, bl_hi, bl_lo
      character(len=STRING_LEN) :: filename_wav, filename_scoeff

      ! Open file.
      fileid = 11
      open(unit=fileid, file=trim(filename), status='new', action='write', &
           form='formatted')

      ! Write header.
      write(fileid, '(a,a)') 'function [J, N, B, bl_scoeff, alpha, scoeff, wcoeff] = ', &
           filename(1:len(trim(filename))-2)
      write(fileid, '(a,a)') '% function [J, N, B, bl_scoeff, alpha, scoeff, wcoeff] = ', &
           filename(1:len(trim(filename))-2)
      write(fileid, '(a)') '%'
      write(fileid, '(a)') '%  Scale discretised wavelet coefficients created by SSHT code'
      write(fileid, '(a)') '%  SSHT code written by Jason McEwen (mcewen@mrao.cam.ac.uk)'
      write(fileid, '(a)') '%'
      write(fileid, '(a)') '%  This file written by: Fortran'
      if(present(comment)) then
         write(fileid, '(a,a)') '%  Comment: ', trim(comment)
      else
         write(fileid, '(a)') '%  Comment: '
      end if
      write(fileid, '(a)') '%'
      write(fileid, '(a)') '%  Variables:'
      write(fileid, '(a)') '%   J: Maximum analysis scale depth.'
      write(fileid, '(a)') '%   B: Harmonic band limit.'
      write(fileid, '(a)') '%   N: Azimuthal band limit.'
      write(fileid, '(a)') '%   bl_scoeff: Upper band limit for scaling coefficients.'
      write(fileid, '(a)') '%   alpha: Basis dilation factor.'
      write(fileid, '(a)') '%   scoeff(1:bl_scoeff,1:bl_scoeff): Scaling coefficients.'
      write(fileid, '(a)') '%   wcoeff{1:J+1}(aa,bb,gg): Wavelet coefficients (ranges of aa, bb and gg, depend on scale).'
      write(fileid, '(a)') 

      ! Write parameters.
      write(fileid, '(a)') '% Size parameters'
      write(fileid, '(a,i27,a)') 'J         = ', J, ';'
      write(fileid, '(a,i27,a)') 'N         = ', N, ';'
      write(fileid, '(a,i27,a)') 'B         = ', B, ';'
      write(fileid, '(a,i27,a)') 'bl_scoeff = ', bl_scoeff, ';'
      write(fileid, '(a,e27.20,a)') 'alpha     = ', alpha, ';'
      write(fileid, '(a)') 

      ! Write wavelet sizes.
      ! Could compute these using bl_lo and bl_hi routines but to avoid any
      ! numerical errors (which anyway shouldn't occur) we get sizes from the wavelet
      ! coefficient arrays.
      write(fileid, '(a)') '% Wavelet coefficient array sizes'
      do jj = 0,J

         ! Set band limits.
         if(ssht_core_assimilate_int(B / (alpha**(jj-1)), bl_hi)) then
            bl_hi = min(bl_hi , B)
         else
            bl_hi = min(ceiling(B / (alpha**(jj-1)) ) , B)
         end if
         if(ssht_core_assimilate_int(B / (alpha**(jj+1)), bl_lo)) then
            bl_lo = max(bl_lo, 0)
         else
            bl_lo = max(floor(B / (alpha**(jj+1)) ), 0)
         end if
         
         ! Write coefficients.
         write(fileid,'(a,i5,a,i15,a)') 'wcoeff_size_aa(', jj+1, ') = ', &
              2*bl_hi-1, ';'
         write(fileid,'(a,i5,a,i15,a)') 'wcoeff_size_bb(', jj+1, ') = ', &
              2*bl_hi, ';'
         write(fileid,'(a,i5,a,i15,a)') 'wcoeff_size_gg(', jj+1, ') = ', &
              N, ';'

      end do
      write(fileid, '(a)') 

      ! Write scaling coefficients to data file.
      write(filename_scoeff, '(a,a)') filename(1:len(trim(filename))-2), '_scoeff.dat'
      fileid_scoeff = 12
      open(unit=fileid_scoeff, file=filename_scoeff, status='new', action='write', &
           form='formatted')

      do el = 0, bl_scoeff-1
         do m = 0, el				
            write(fileid_scoeff, '(e27.20)') real(scoeff(el, m),dp)
         end do
      end do
      do el = 0, bl_scoeff-1
         do m = 0, el				
            write(fileid_scoeff, '(e27.20)') real(aimag(scoeff(el, m)),dp)
         end do
      end do

      ! Write matlab code to read scaling coefficients.
      write(fileid, '(a)') '% Scaling coefficients'
      write(fileid, '(a,a,a)') 'scoeff_array = load(''', trim(filename_scoeff), ''');'
      write(fileid, '(a)') 'n = length(scoeff_array) / 2;'
      write(fileid, '(a)') 'scoeff_array = scoeff_array(1:n) + i*scoeff_array(n+1:end);'
      write(fileid, '(a)') 'scoeff = zeros(bl_scoeff,bl_scoeff);'
      write(fileid, '(a)') 'ind_start = 1;'
      write(fileid, '(a)') 'for el = 0:bl_scoeff-1'
      write(fileid, '(a)') '    ind_end = ind_start - 1 + el + 1;'
      write(fileid, '(a)') '    scoeff(el+1, 1:el+1) = scoeff_array(ind_start:ind_end).'';'
      write(fileid, '(a)') '    ind_start = ind_end + 1;'
      write(fileid, '(a)') 'end'
      write(fileid, '(a)') 

      ! Write wavelet coefficients to data file.
      write(filename_wav, '(a,a)') filename(1:len(trim(filename))-2), '_wcoeff.dat'
      fileid_wav = 13
      open(unit=fileid_wav, file=filename_wav, status='new', action='write', &
           form='formatted')

      do jj = 0,J    

         ! Set band limits.
         if(ssht_core_assimilate_int(B / (alpha**(jj-1)), bl_hi)) then
            bl_hi = min(bl_hi , B)
         else
            bl_hi = min(ceiling(B / (alpha**(jj-1)) ) , B)
         end if
         if(ssht_core_assimilate_int(B / (alpha**(jj+1)), bl_lo)) then
            bl_lo = max(bl_lo, 0)
         else
            bl_lo = max(floor(B / (alpha**(jj+1)) ), 0)
         end if

         ! Write coefficients.
         do gg = 0, N-1
            do bb = 0, 2*bl_hi-1
               do aa = 0, 2*bl_hi-2
                  write(fileid_wav,'(e27.20)') wav(jj,aa,bb,gg)
               end do
            end do
         end do

      end do

      close(fileid_wav)
   
      ! Write matlab code to read wavelet coefficients.
      write(fileid, '(a)') '% Wavelet coefficients'
      write(fileid, '(a,a,a)') 'wcoeff_array = load(''', trim(filename_wav), ''');'
      write(fileid, '(a)') 'ind_start = 1;'
      write(fileid, '(a)') 'for jj = 1:J+1'
      write(fileid, '(a)') '    ind_end = ind_start - 1 + wcoeff_size_aa(jj) * wcoeff_size_bb(jj) * wcoeff_size_gg(jj);'
      write(fileid, '(a,a)') '    wcoeff{jj} = reshape(wcoeff_array(ind_start:ind_end),', &
           'wcoeff_size_aa(jj), wcoeff_size_bb(jj), wcoeff_size_gg(jj));'
      write(fileid, '(a)') '    ind_start = ind_end+1;'
      write(fileid, '(a)') 'end'
      write(fileid, '(a)') 

      ! Close file.
      close(fileid)

    end subroutine ssht_fileio_matlab_wav_write_static


    !--------------------------------------------------------------------------
    ! ssht_fileio_matlab_wav_write_dynamic
    !
    !! Writes (dynamically allocated) wavelet and scaling coefficients to an 
    !! output SSHT formatted .m matlab file ad corresponding data .dat file.
    !!
    !! Variables:
    !!  - wavdyn(0:J)%coeff: Dynamically allocated wavelet coefficients for
    !!    each scale [input].
		!!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
		!!    [input].
		!!  - J: Maximum analysis scale depth [input].
		!!  - B: Harmonic band limit [input].
		!!  - N: Azimuthal band limit [input].
		!!  - bl_scoeff: Upper band limit for scaling coefficients [input].
		!!  - alpha: Basis dilation factor [input].
    !!  - filename: Name of the output matlab file to write [input].
    !!  - [comment]: Optional comment string to be added to the output matlab 
    !!    file header if present [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - June 2008
    !
    ! Revisions:
    !   June 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_matlab_wav_write_dynamic(wavdyn, scoeff, J, B, N, &
      bl_scoeff, alpha, filename, comment)

			integer, intent(in) :: J
			integer, intent(in) :: B
			integer, intent(in) :: N
			integer, intent(in) :: bl_scoeff
			real(dp), intent(in) :: alpha
      type(ssht_wav_abg), intent(in), allocatable :: wavdyn(:)
			complex(dpc), intent(in) :: scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      integer :: fileid, fileid_wav, fileid_scoeff
      integer :: el, m, jj, aa, bb, gg
      character(len=STRING_LEN) :: filename_wav, filename_scoeff

      ! Open file.
      fileid = 11
      open(unit=fileid, file=trim(filename), status='new', action='write', &
           form='formatted')

      ! Write header.
      write(fileid, '(a,a)') 'function [J, N, B, bl_scoeff, alpha, scoeff, wcoeff] = ', &
           filename(1:len(trim(filename))-2)
      write(fileid, '(a,a)') '% function [J, N, B, bl_scoeff, alpha, scoeff, wcoeff] = ', &
           filename(1:len(trim(filename))-2)
      write(fileid, '(a)') '%'
      write(fileid, '(a)') '%  Scale discretised wavelet coefficients created by SSHT code'
      write(fileid, '(a)') '%  SSHT code written by Jason McEwen (mcewen@mrao.cam.ac.uk)'
      write(fileid, '(a)') '%'
      write(fileid, '(a)') '%  This file written by: Fortran'
      if(present(comment)) then
         write(fileid, '(a,a)') '%  Comment: ', trim(comment)
      else
         write(fileid, '(a)') '%  Comment: '
      end if
      write(fileid, '(a)') '%'
      write(fileid, '(a)') '%  Variables:'
      write(fileid, '(a)') '%   J: Maximum analysis scale depth.'
      write(fileid, '(a)') '%   B: Harmonic band limit.'
      write(fileid, '(a)') '%   N: Azimuthal band limit.'
      write(fileid, '(a)') '%   bl_scoeff: Upper band limit for scaling coefficients.'
      write(fileid, '(a)') '%   alpha: Basis dilation factor.'
      write(fileid, '(a)') '%   scoeff(1:bl_scoeff,1:bl_scoeff): Scaling coefficients.'
      write(fileid, '(a)') '%   wcoeff{1:J+1}(aa,bb,gg): Wavelet coefficients (ranges of aa, bb and gg, depend on scale).'
      write(fileid, '(a)') 

      ! Write parameters.
      write(fileid, '(a)') '% Size parameters'
      write(fileid, '(a,i27,a)') 'J         = ', J, ';'
      write(fileid, '(a,i27,a)') 'N         = ', N, ';'
      write(fileid, '(a,i27,a)') 'B         = ', B, ';'
      write(fileid, '(a,i27,a)') 'bl_scoeff = ', bl_scoeff, ';'
      write(fileid, '(a,e27.20,a)') 'alpha     = ', alpha, ';'
      write(fileid, '(a)') 

      ! Write wavelet sizes.
      ! Could compute these using bl_lo and bl_hi routines but to avoid any
      ! numerical errors (which anyway shouldn't occur) we get sizes from the wavelet
      ! coefficient arrays.
      write(fileid, '(a)') '% Wavelet coefficient array sizes'
      do jj = 0,J
         write(fileid,'(a,i5,a,i15,a)') 'wcoeff_size_aa(', jj+1, ') = ', &
              size(wavdyn(jj)%coeff,1), ';'
         write(fileid,'(a,i5,a,i15,a)') 'wcoeff_size_bb(', jj+1, ') = ', &
              size(wavdyn(jj)%coeff,2), ';'
         write(fileid,'(a,i5,a,i15,a)') 'wcoeff_size_gg(', jj+1, ') = ', &
              size(wavdyn(jj)%coeff,3), ';'
      end do
      write(fileid, '(a)') 

      ! Write scaling coefficients to data file.
      write(filename_scoeff, '(a,a)') filename(1:len(trim(filename))-2), '_scoeff.dat'
      fileid_scoeff = 12
      open(unit=fileid_scoeff, file=filename_scoeff, status='new', action='write', &
           form='formatted')

      do el = 0, bl_scoeff-1
         do m = 0, el				
            write(fileid_scoeff, '(e27.20)') real(scoeff(el, m),dp)
         end do
      end do
      do el = 0, bl_scoeff-1
         do m = 0, el				
            write(fileid_scoeff, '(e27.20)') real(aimag(scoeff(el, m)),dp)
         end do
      end do

      close(fileid_scoeff)

      ! Write matlab code to read scaling coefficients.
      write(fileid, '(a)') '% Scaling coefficients'
      write(fileid, '(a,a,a)') 'scoeff_array = load(''', trim(filename_scoeff), ''');'
      write(fileid, '(a)') 'n = length(scoeff_array) / 2;'
      write(fileid, '(a)') 'scoeff_array = scoeff_array(1:n) + i*scoeff_array(n+1:end);'
      write(fileid, '(a)') 'scoeff = zeros(bl_scoeff,bl_scoeff);'
      write(fileid, '(a)') 'ind_start = 1;'
      write(fileid, '(a)') 'for el = 0:bl_scoeff-1'
      write(fileid, '(a)') '    ind_end = ind_start - 1 + el + 1;'
      write(fileid, '(a)') '    scoeff(el+1, 1:el+1) = scoeff_array(ind_start:ind_end).'';'
      write(fileid, '(a)') '    ind_start = ind_end + 1;'
      write(fileid, '(a)') 'end'
      write(fileid, '(a)') 

      ! Write wavelet coefficients to data file.
      write(filename_wav, '(a,a)') filename(1:len(trim(filename))-2), '_wcoeff.dat'
      fileid_wav = 13
      open(unit=fileid_wav, file=filename_wav, status='new', action='write', &
           form='formatted')

      do jj = 0,J    

         do gg = 0, size(wavdyn(jj)%coeff,3)-1
            do bb = 0, size(wavdyn(jj)%coeff,2)-1
               do aa = 0, size(wavdyn(jj)%coeff,1)-1
                  write(fileid_wav,'(e27.20)') wavdyn(jj)%coeff(aa,bb,gg)
               end do
            end do
         end do

      end do

      close(fileid_wav)
   
      ! Write matlab code to read wavelet coefficients.
      write(fileid, '(a)') '% Wavelet coefficients'
      write(fileid, '(a,a,a)') 'wcoeff_array = load(''', trim(filename_wav), ''');'
      write(fileid, '(a)') 'ind_start = 1;'
      write(fileid, '(a)') 'for jj = 1:J+1'
      write(fileid, '(a)') '    ind_end = ind_start - 1 + wcoeff_size_aa(jj) * wcoeff_size_bb(jj) * wcoeff_size_gg(jj);'
      write(fileid, '(a,a)') '    wcoeff{jj} = reshape(wcoeff_array(ind_start:ind_end),', &
           'wcoeff_size_aa(jj), wcoeff_size_bb(jj), wcoeff_size_gg(jj));'
      write(fileid, '(a)') '    ind_start = ind_end+1;'
      write(fileid, '(a)') 'end'
      write(fileid, '(a)') 

      ! Close file.
      close(fileid)

    end subroutine ssht_fileio_matlab_wav_write_dynamic


    !--------------------------------------------------------------------------
    ! ssht_fileio_matlab_wav_read_static
    !
    !! Reads (statically allocated) wavelet and scaling coefficients from a
    !! SSHT formatted .m matlab file and corresponding .dat data files.
    !!
    !! Notes:
    !!   - Memory for the wavelet and scaling coefficients is allocated herein
		!!     and should be freed by the calling routine.
    !!   - Input file must be precisely formatted as written by appropriate 
    !!     SSHT Fortran or Matlab routine (file error checking is *not* 
    !!     performed).
    !!
    !! Variables:
    !!  - wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1): Wavelet coefficients [output].
		!!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
		!!    [output].
		!!  - J: Maximum analysis scale depth [output].
		!!  - B: Harmonic band limit [output].
		!!  - N: Azimuthal band limit [output].
		!!  - bl_scoeff: Upper band limit for scaling coefficients [output].
		!!  - alpha: Basis dilation factor [output].
    !!  - filename: Name of the matlab input file containing the wavelet
    !!    and scaling coefficients to be read [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - May 2008
    !
    ! Revisions:
    !   May 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_matlab_wav_read_static(wav, scoeff, J, B, N, &
         bl_scoeff, alpha, filename)

			integer, intent(out) :: J
			integer, intent(out) :: B
			integer, intent(out) :: N
			integer, intent(out) :: bl_scoeff
			real(dp), intent(out) :: alpha
			real(dp), allocatable, intent(out) :: wav(:,:,:,:)
			complex(dpc), allocatable, intent(out) :: scoeff(:,:)
      character(len=*), intent(in) :: filename

      integer :: fileid, fileid_wav, fileid_scoeff
      integer :: el, m, jj, aa, bb, gg
      character(len=STRING_LEN) :: line, line_trun
      integer :: iline, fail = 0
      integer, allocatable :: size1(:), size2(:), size3(:)
      real(dp) :: scoeff_real, scoeff_imag
      character(len=STRING_LEN) :: filename_wav, filename_scoeff

      ! Open file.
      fileid = 12
      open(unit=fileid, file=trim(filename), status='old', action='read', &
           form='formatted')

      ! Skip first 19 lines.
      do iline = 1,19
         read(fileid, '(a)') line
      end do

      ! Read parameters
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) J
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) N
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) B
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) bl_scoeff
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) alpha

      ! Skip 2 lines.
      do iline = 1,2
         read(fileid, '(a)') line
      end do

      ! Allocate space for wavelet and scaling coefficients.
      fail = 0
			allocate(wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1), stat=fail)
			allocate(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
			if(fail /= 0) then
				call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_fileio_matlab_wav_read_static')
			end if
      wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1) = 0d0
			scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = 0d0

      ! Read wavelet coefficient sizes.
      allocate(size1(0:J), stat=fail)
      allocate(size2(0:J), stat=fail)
      allocate(size3(0:J), stat=fail)
			if(fail /= 0) then
				call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_fileio_matlab_wav_read_static')
			end if
      do jj = 0,J

         read(fileid, '(a)') line
         line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
         read(line_trun,*) size1(jj)
         read(fileid, '(a)') line
         line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
         read(line_trun,*) size2(jj)
         read(fileid, '(a)') line
         line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
         read(line_trun,*) size3(jj)

      end do

      close(fileid)

      ! Read scaling coefficients.
      write(filename_scoeff, '(a,a)') filename(1:len(trim(filename))-2), '_scoeff.dat'
      fileid_scoeff = 13
      open(unit=fileid_scoeff, file=filename_scoeff, status='old', action='read', &
           form='formatted')

      do el = 0, bl_scoeff-1
         do m = 0, el			
            read(fileid_scoeff, '(e27.20)') scoeff_real            
            scoeff(el,m) = scoeff_real
         end do
      end do
      do el = 0, bl_scoeff-1
         do m = 0, el				
            read(fileid_scoeff, '(e27.20)') scoeff_imag
            scoeff(el,m) = cmplx(real(scoeff(el,m),dp), scoeff_imag)
         end do
      end do

      close(fileid_scoeff)

      ! Read wavelet coefficients.
      write(filename_wav, '(a,a)') filename(1:len(trim(filename))-2), '_wcoeff.dat'
      fileid_wav = 14
      open(unit=fileid_wav, file=filename_wav, status='old', action='read', &
           form='formatted')

      do jj = 0,J     
         do gg = 0, size3(jj)-1
            do bb = 0, size2(jj)-1
               do aa = 0, size1(jj)-1
                  read(fileid_wav, '(e27.20)') wav(jj,aa,bb,gg)     
               end do
            end do
         end do

      end do

      close(fileid_wav)

      ! Free memory.
      deallocate(size1, size2, size3)

    end subroutine ssht_fileio_matlab_wav_read_static


    !--------------------------------------------------------------------------
    ! ssht_fileio_matlab_wav_read_dynamic
    !
    !! Reads (dynamically allocated) wavelet and scaling coefficients from a
    !! SSHT formatted .m matlab file and corresponding .dat data files.
    !!
    !! Notes:
    !!   - Memory for the wavelet and scaling coefficients is allocated herein
		!!     and should be freed by the calling routine.
    !!   - Input file must be precisely formatted as written by appropriate 
    !!     SSHT Fortran or Matlab routine (file error checking is *not* 
    !!     performed).
    !!
    !! Variables:
    !!  - wavdyn(0:J)%coeff: Dynamically allocated wavelet coefficients for
    !!    each scale (memory allocated herein) [output].
		!!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
		!!    [output].
		!!  - J: Maximum analysis scale depth [output].
		!!  - B: Harmonic band limit [output].
		!!  - N: Azimuthal band limit [output].
		!!  - bl_scoeff: Upper band limit for scaling coefficients [output].
		!!  - alpha: Basis dilation factor [output].
    !!  - filename: Name of the matlab input file containing the wavelet
    !!    and scaling coefficients to be read [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - May 2008
    !
    ! Revisions:
    !   May 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_matlab_wav_read_dynamic(wavdyn, scoeff, J, B, N, &
         bl_scoeff, alpha, filename)

			integer, intent(out) :: J
			integer, intent(out) :: B
			integer, intent(out) :: N
			integer, intent(out) :: bl_scoeff
			real(dp), intent(out) :: alpha
      type(ssht_wav_abg), intent(out), allocatable :: wavdyn(:)
			complex(dpc), allocatable, intent(out) :: scoeff(:,:)
      character(len=*), intent(in) :: filename

      integer :: fileid, fileid_wav, fileid_scoeff
      integer :: el, m, jj, aa, bb, gg
      character(len=STRING_LEN) :: line, line_trun
      integer :: iline, fail = 0
      integer :: size1, size2, size3
      real(dp) :: scoeff_real, scoeff_imag
      character(len=STRING_LEN) :: filename_wav, filename_scoeff

      ! Open file.
      fileid = 12
      open(unit=fileid, file=trim(filename), status='old', action='read', &
           form='formatted')

      ! Skip first 19 lines.
      do iline = 1,19
         read(fileid, '(a)') line
      end do

      ! Read parameters
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) J
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) N
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) B
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) bl_scoeff
      read(fileid, '(a)') line
      line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
      read(line_trun,*) alpha

      ! Skip 2 lines.
      do iline = 1,2
         read(fileid, '(a)') line
      end do

      ! Allocate space for wavelet and scaling coefficients
      ! (while reading coefficient array sizes).
      fail = 0
      allocate(wavdyn(0:J), stat=fail)
			allocate(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
			if(fail /= 0) then
				call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_fileio_matlab_wav_read_dynamic')
			end if
			scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = 0d0

      do jj = 0,J

         read(fileid, '(a)') line
         line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
         read(line_trun,*) size1
         read(fileid, '(a)') line
         line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
         read(line_trun,*) size2
         read(fileid, '(a)') line
         line_trun = trim(line(index(line,'=')+1:len(trim(line))-1))
         read(line_trun,*) size3

         ! Allocate space for wavelet coefficients at given scale.
        allocate(wavdyn(jj)%coeff(0:size1-1, 0:size2-1, 0:size3-1), stat=fail)
        if(fail /= 0) then
           call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, &
                'ssht_fileio_matlab_wav_read_dynamic')
        end if
        wavdyn(jj)%coeff(0:size1-1, 0:size2-1, 0:size3-1) = 0d0

      end do

      close(fileid)

      ! Read scaling coefficients.
      write(filename_scoeff, '(a,a)') filename(1:len(trim(filename))-2), '_scoeff.dat'
      fileid_scoeff = 13
      open(unit=fileid_scoeff, file=filename_scoeff, status='old', action='read', &
           form='formatted')

      do el = 0, bl_scoeff-1
         do m = 0, el			
            read(fileid_scoeff, '(e27.20)') scoeff_real            
            scoeff(el,m) = scoeff_real
         end do
      end do
      do el = 0, bl_scoeff-1
         do m = 0, el				
            read(fileid_scoeff, '(e27.20)') scoeff_imag
            scoeff(el,m) = cmplx(real(scoeff(el,m),dp), scoeff_imag)
         end do
      end do

      close(fileid_scoeff)

      ! Read wavelet coefficients.
      write(filename_wav, '(a,a)') filename(1:len(trim(filename))-2), '_wcoeff.dat'
      fileid_wav = 14
      open(unit=fileid_wav, file=filename_wav, status='old', action='read', &
           form='formatted')

      do jj = 0,J     
         do gg = 0, size(wavdyn(jj)%coeff,3)-1
            do bb = 0, size(wavdyn(jj)%coeff,2)-1
               do aa = 0, size(wavdyn(jj)%coeff,1)-1
                  read(fileid_wav, '(e27.20)') wavdyn(jj)%coeff(aa,bb,gg)       
               end do
            end do
         end do

      end do

      close(fileid_wav)

    end subroutine ssht_fileio_matlab_wav_read_dynamic


    !--------------------------------------------------------------------------
    ! Fits file IO
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    ! ssht_fileio_fits_wav_write_static
    !
    !! Writes (statically allocated) wavelet and scaling coefficients to an
    !! output SSHT formatted fits file. 
    !!
    !! Variables:
    !!  - wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1): Wavelet coefficients [input].
		!!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
		!!    [input].
		!!  - J: Maximum analysis scale depth [input].
		!!  - B: Harmonic band limit [input].
		!!  - N: Azimuthal band limit [input].
		!!  - bl_scoeff: Upper band limit for scaling coefficients [input].
		!!  - alpha: Basis dilation factor [input].
    !!  - filename: Name of the output fits file to write [input].
    !!  - [comment]: Optional comment string to be added to the output fits 
    !!    file header if present [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2007
    !
    ! Revisions:
    !   November 2007 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_fits_wav_write_static(wav, scoeff, J, B, N, &
         bl_scoeff, alpha, filename, comment)

			integer, intent(in) :: J
			integer, intent(in) :: B
			integer, intent(in) :: N
			integer, intent(in) :: bl_scoeff
			real(dp), intent(in) :: alpha
			real(dp), intent(in) :: wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1)
			complex(dpc), intent(in) :: scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      integer :: status,unit,blocksize,bitpix
      integer :: group,dim1,dim2
      logical :: simple, extend, file_exists
      integer :: decimals
      integer :: naxis
      integer :: naxes(1)
			integer :: naxes_wcoeff(1:3)
			integer :: naxes_scoeff(1:2)
      integer :: jj, bl_hi, bl_lo

      ! Define FITS parameters.
      bitpix=-64    ! Real double precision.
      status=0     ! Initialse error status to zero.
      decimals = 8 ! Number of decimals used to store alpha.

      ! Check if file already exists.
      call ssht_fileio_fits_exists(filename, status, file_exists)
      if(file_exists) then
         call ssht_error(SSHT_ERROR_FILEIO, &
						'ssht_fileio_fits_wav_write', &
						comment_add='File already exists')
        stop
      end if

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Create the new empty fits file.
      blocksize=1  ! Historical artifact that is ignored.
      call ftinit(unit,filename,blocksize,status)

      ! Write primary header.
      simple=.true.
      extend=.true.
      naxis=0
      naxes(1)=0
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      ! Write additional header keywords.
      call ftpcom(unit, &
        '  Scale discretised wavelet coefficients created by SSHT code',status)
      call ftpcom(unit, &
        '  Written by Jason McEwen (mcewen@mrao.cam.ac.uk)',status)
      call ftpcom(unit, &
        '  Primary extension empty.',status)
      call ftpdat(unit,status)    ! Add date
      if(present(comment)) then 
         call ftpcom(unit, comment, status)
      end if
      call ftpkyj(unit,'J', J, 'No. of scales',status)
			call ftpkyj(unit,'B', B, 'Band limit',status)
			call ftpkyj(unit,'N', N, 'Azimuthal band limit',status)
			call ftpkyj(unit,'BLSCOEFF', bl_scoeff, 'Scaling coefficient band limit',status)
			call ftpkyd(unit,'ALPHA', alpha,decimals, 'Alpha',status)

			! Write scaling coefficients to image extensions.			
			naxis=2
			naxes_scoeff(1) = bl_scoeff
			naxes_scoeff(2) = bl_scoeff
			group = 1
			dim1 = bl_scoeff
			! Write real part of scaling coefficients.
			call ftiimg(unit,bitpix,naxis,naxes_scoeff,status)
			call ftpkys(unit,'EXTNAME','SCOEFF_REAL','entension name',status)
			call ftp2dd(unit, group, dim1, bl_scoeff, bl_scoeff, &
				real(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1),dp), status)
			! Write imaginary part of scaling coefficients.
			call ftiimg(unit,bitpix,naxis,naxes_scoeff,status)
			call ftpkys(unit,'EXTNAME','SCOEFF_IMAG','entension name',status)
			call ftp2dd(unit, group, dim1, bl_scoeff, bl_scoeff, &
				real(aimag(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)),dp), status)

      ! Insert additional image extensions for wavelet coefficients at 
      ! each scale.
      do jj = 0,J

        ! Set band limits.
				if(ssht_core_assimilate_int(B / (alpha**(jj-1)), bl_hi)) then
					bl_hi = min(bl_hi , B)
				else
					bl_hi = min(ceiling(B / (alpha**(jj-1)) ) , B)
				end if
				if(ssht_core_assimilate_int(B / (alpha**(jj+1)), bl_lo)) then
					bl_lo = max(bl_lo, 0)
				else
					bl_lo = max(floor(B / (alpha**(jj+1)) ), 0)
				end if

	      naxis=3
				naxes_wcoeff(1) = 2*bl_hi-1
				naxes_wcoeff(2) = 2*bl_hi
				naxes_wcoeff(3) = N

        ! Insert a new image extension.
        call ftiimg(unit,bitpix,naxis,naxes_wcoeff,status)

        ! Write additional header keywords.
        call ftpkys(unit,'EXTNAME','WAV','entension name',status)

        ! Write wavelet coefficients for particular scale as a 3D data cube.
        group = 1
        dim1 = 2*bl_hi-1
        dim2 = 2*bl_hi
        call ftp3dd(unit, group, dim1, dim2, 2*bl_hi-1, 2*bl_hi, N, &
					wav(jj, 0:2*bl_hi-2, 0:2*bl_hi-1, 0:N-1), status)

      end do

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call ssht_fileio_fits_error_check(status, .true.)

    end subroutine ssht_fileio_fits_wav_write_static


    !--------------------------------------------------------------------------
    ! ssht_fileio_fits_wav_write_dynamic
    !
    !! Writes (dynamically allocated) wavelet and scaling coefficients to an output SSHT formatted 
		!! fits file. 
    !!
    !! Variables:
    !!  - wavdyn(0:J)%coeff: Dynamically allocated wavelet coefficients for
    !!    each scale [input].
		!!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
		!!    [input].
		!!  - J: Maximum analysis scale depth [input].
		!!  - B: Harmonic band limit [input].
		!!  - N: Azimuthal band limit [input].
		!!  - bl_scoeff: Upper band limit for scaling coefficients [input].
		!!  - alpha: Basis dilation factor [input].
    !!  - filename: Name of the output fits file to write [input].
    !!  - [comment]: Optional comment string to be added to the output fits 
    !!    file header if present [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - February 2008
    !
    ! Revisions:
    !   February 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_fits_wav_write_dynamic(wavdyn, scoeff, J, B, N, &
      bl_scoeff, alpha, filename, comment)

			integer, intent(in) :: J
			integer, intent(in) :: B
			integer, intent(in) :: N
			integer, intent(in) :: bl_scoeff
			real(dp), intent(in) :: alpha
      type(ssht_wav_abg), intent(in), allocatable :: wavdyn(:)
			complex(dpc), intent(in) :: scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)
      character(len=*), intent(in) :: filename
      character(len=*), intent(in), optional :: comment

      integer :: status,unit,blocksize,bitpix
      integer :: group,dim1,dim2
      logical :: simple, extend, file_exists
      integer :: decimals
      integer :: naxis
      integer :: naxes(1)
			integer :: naxes_wcoeff(1:3)
			integer :: naxes_scoeff(1:2)
      integer :: jj, bl_hi, bl_lo

      ! Define FITS parameters.
      bitpix=-64    ! Real double precision.
      status=0     ! Initialse error status to zero.
      decimals = 8 ! Number of decimals used to store alpha.

      ! Check if file already exists.
      call ssht_fileio_fits_exists(filename, status, file_exists)
      if(file_exists) then
         call ssht_error(SSHT_ERROR_FILEIO, &
						'ssht_fileio_fits_wav_write_dynamic', &
						comment_add='File already exists')
        stop
      end if

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Create the new empty fits file.
      blocksize=1  ! Historical artifact that is ignored.
      call ftinit(unit,filename,blocksize,status)

      ! Write primary header.
      simple=.true.
      extend=.true.
      naxis=0
      naxes(1)=0
      call ftphpr(unit,simple,bitpix,naxis,naxes,0,1,extend,status)

      ! Write additional header keywords.
      call ftpcom(unit, &
        '  Scale discretised wavelet coefficients created by SSHT code',status)
      call ftpcom(unit, &
        '  Written by Jason McEwen (mcewen@mrao.cam.ac.uk)',status)
      call ftpcom(unit, &
        '  Primary extension empty.',status)
      call ftpdat(unit,status)    ! Add date
      if(present(comment)) then 
         call ftpcom(unit, comment, status)
      end if
      call ftpkyj(unit,'J', J, 'No. of scales',status)
			call ftpkyj(unit,'B', B, 'Band limit',status)
			call ftpkyj(unit,'N', N, 'Azimuthal band limit',status)
			call ftpkyj(unit,'BLSCOEFF', bl_scoeff, 'Scaling coefficient band limit',status)
			call ftpkyd(unit,'ALPHA', alpha,decimals, 'Alpha',status)

			! Write scaling coefficients to image extensions.			
			naxis=2
			naxes_scoeff(1) = bl_scoeff
			naxes_scoeff(2) = bl_scoeff
			group = 1
			dim1 = bl_scoeff
			! Write real part of scaling coefficients.
			call ftiimg(unit,bitpix,naxis,naxes_scoeff,status)
			call ftpkys(unit,'EXTNAME','SCOEFF_REAL','entension name',status)
			call ftp2dd(unit, group, dim1, bl_scoeff, bl_scoeff, &
				real(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1),dp), status)
			! Write imaginary part of scaling coefficients.
			call ftiimg(unit,bitpix,naxis,naxes_scoeff,status)
			call ftpkys(unit,'EXTNAME','SCOEFF_IMAG','entension name',status)
			call ftp2dd(unit, group, dim1, bl_scoeff, bl_scoeff, &
				real(aimag(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1)),dp), status)

      ! Insert additional image extensions for wavelet coefficients at 
      ! each scale.
      do jj = 0,J

        naxis=3
        naxes_wcoeff(1) = size(wavdyn(jj)%coeff, 1)
        naxes_wcoeff(2) = size(wavdyn(jj)%coeff, 2)
        naxes_wcoeff(3) = size(wavdyn(jj)%coeff, 3)

        ! Insert a new image extension.
        call ftiimg(unit,bitpix,naxis,naxes_wcoeff,status)

        ! Write additional header keywords.
        call ftpkys(unit,'EXTNAME','WAV','entension name',status)

        ! Write wavelet coefficients for particular scale as a 3D data cube.
        group = 1
        dim1 = size(wavdyn(jj)%coeff, 1)
        dim2 = size(wavdyn(jj)%coeff, 2)
        call ftp3dd(unit, group, dim1, dim2, &
             naxes_wcoeff(1), naxes_wcoeff(2), naxes_wcoeff(3), &
             wavdyn(jj)%coeff(0:size(wavdyn(jj)%coeff,1)-1, &
                        0:size(wavdyn(jj)%coeff,2)-1, &
                        0:size(wavdyn(jj)%coeff,3)-1), status)

      end do

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call ssht_fileio_fits_error_check(status, .true.)

    end subroutine ssht_fileio_fits_wav_write_dynamic


    !--------------------------------------------------------------------------
    ! ssht_fileio_fits_wav_read_static
    !
    !! Reads (statically allocated) wavelet and scaling coefficients from a
    !! SSHT formatted fits file.
    !!
    !! Notes:
    !!   - Memory for the wavelet and scaling coefficients is allocated herein
		!!     and should be freed by the calling routine.
    !!
    !! Variables:
    !!  - wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1): Wavelet coefficients [output].
		!!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
		!!    [output].
		!!  - J: Maximum analysis scale depth [output].
		!!  - B: Harmonic band limit [output].
		!!  - N: Azimuthal band limit [output].
		!!  - bl_scoeff: Upper band limit for scaling coefficients [output].
		!!  - alpha: Basis dilation factor [output].
    !!  - filename: Name of the fits input file containing the wavelet
    !!    and scaling coefficients to be read [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2007
    !
    ! Revisions:
    !   November 2007 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_fits_wav_read_static(wav, scoeff, J, B, N, &
         bl_scoeff, alpha, filename)

			integer, intent(out) :: J
			integer, intent(out) :: B
			integer, intent(out) :: N
			integer, intent(out) :: bl_scoeff
			real(dp), intent(out) :: alpha
			real(dp), allocatable, intent(out) :: wav(:,:,:,:)
			complex(dpc), allocatable, intent(out) :: scoeff(:,:)
      character(len=*), intent(in) :: filename

      character(len=20) :: comment
      integer :: status, unit, blocksize, readwrite
      integer :: nhdu, hdutype, naxis1, naxis2, naxis3
      integer :: hdunum, hdunum_check
      logical :: anynull, file_exists
      integer :: group, dim1, dim2
      real(sp) :: nullval
			integer :: jj, bl_lo, bl_hi, fail
			real(dp), allocatable :: scoeff_temp(:,:)

			! Initialse error status to zero.
      status = 0   

      ! Check if file already exists.
      call ssht_fileio_fits_exists(filename, status, file_exists)
      if(.not. file_exists) then
         call ssht_error(SSHT_ERROR_FILEIO, &
           'ssht_fileio_fits_wav_read', &
           comment_add='File does not exist')
      end if

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Open file as readonly. 
      readwrite = 0    ! Open as readonly.
      call ftopen(unit, filename, readwrite, blocksize, status)

      ! Read size of image from variable size keywords.
      call ftgkyj(unit, 'J', J, comment, status)
      call ftgkyj(unit, 'B', B, comment, status)
      call ftgkyj(unit, 'N', N, comment, status)
      call ftgkyj(unit, 'BLSCOEFF', bl_scoeff, comment, status)
			 call ftgkyd(unit, 'ALPHA', alpha, comment, status)

      ! Check correct number of HDUs in input file.
      ! (First two extension due to primary and dilation, then have one
      ! extension for coefficients at each scale.)
      hdunum = 1 + 2 + (J+1)
      call ftthdu(unit, hdunum_check, status)  ! Number extensions in file.
      if(hdunum_check /= hdunum) then
         call ssht_error(SSHT_ERROR_FILEIO, &
           'ssht_fileio_fits_wav_read', &
           comment_add='Invalid number of headers')
      end if

			! Allocate memory.
      fail = 0
			allocate(wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1), stat=fail)
			allocate(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
			allocate(scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
			if(fail /= 0) then
				call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_fileio_fits_wav_read')
			end if
			wav(0:J, 0:2*B-2, 0:2*B-1, 0:N-1) = 0d0
			scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = 0d0

			! Read real part of scaling coefficients.
			nhdu = 2 
			call ftmahd(unit, nhdu, hdutype, status)  
			call ftgkyj(unit, 'NAXIS1', naxis1, comment, status)
			call ftgkyj(unit, 'NAXIS2', naxis2, comment, status)
			if(naxis1 /= bl_scoeff .or. naxis2 /= bl_scoeff) then
				call ssht_error(SSHT_ERROR_FILEIO, &
					'ssht_fileio_fits_wav_read', &
					comment_add='Inconsistent scaling coefficients sizes in file')
			end if
			group = 1
			nullval = -999
			dim1 = bl_scoeff
			call ftg2dd(unit, group, nullval, dim1, bl_scoeff, bl_scoeff, &
				scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1), anynull, status)
			scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1)

			! Read imaginary part of scaling coefficients.
			nhdu = 3  
			call ftmahd(unit, nhdu, hdutype, status)  
			call ftgkyj(unit, 'NAXIS1', naxis1, comment, status)
			call ftgkyj(unit, 'NAXIS2', naxis2, comment, status)
			if(naxis1 /= bl_scoeff .or. naxis2 /= bl_scoeff) then
				call ssht_error(SSHT_ERROR_FILEIO, &
					'ssht_fileio_fits_wav_read', &
					comment_add='Inconsistent scaling coefficients sizes in file')
			end if
			group = 1
			nullval = -999
			dim1 = bl_scoeff
			call ftg2dd(unit, group, nullval, dim1, bl_scoeff, bl_scoeff, &
				scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1), anynull, status)
			scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) &
				+ I * scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1)

      ! Read wavelet coefficients at each scale.
      do jj = 0,J
	
				! Move to next extension.
				nhdu = nhdu + 1 
				call ftmahd(unit, nhdu, hdutype, status)  

				! Read header and check correct sizes.
				call ftgkyj(unit, 'NAXIS1', naxis1, comment, status)
				call ftgkyj(unit, 'NAXIS2', naxis2, comment, status)
				call ftgkyj(unit, 'NAXIS3', naxis3, comment, status)

				! Read coefficients as 3D data cube.
				group = 1
				nullval = -999
				dim1 = naxis1
				dim2 = naxis2
				call ftg3dd(unit, group, nullval, dim1, dim2, &
					naxis1, naxis2, N, &
					wav(jj, 0:naxis1-1, 0:naxis2-1, 0:naxis3-1), anynull, status)
         if(anynull) then
            call ssht_error(SSHT_ERROR_FILEIO, &
              'ssht_fileio_fits_wav_read', &
              comment_add='Null wavelet coefficients read from file')
         end if

      end do

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call ssht_fileio_fits_error_check(status, .true.)

      ! Free temporary storage used.
      deallocate(scoeff_temp)

    end subroutine ssht_fileio_fits_wav_read_static


    !--------------------------------------------------------------------------
    ! ssht_fileio_fits_wav_read_dynamic
    !
    !! Reads (dynamically allocated) wavelet and scaling coefficients from a
    !! SSHT formatted fits file.
    !!
    !! Notes:
    !!   - Memory for the wavelet and scaling coefficients is allocated herein
		!!     and should be freed by the calling routine.
    !!
    !! Variables:
    !!  - wavdyn(0:J)%coeff: Dynamically allocated wavelet coefficients for
    !!    each scale (memory allocated herein) [output].
		!!  - scoeff(0:bl_scoeff-1, 0:bl_scoeff-1): Scaling coefficients 
		!!    [output].
		!!  - J: Maximum analysis scale depth [output].
		!!  - B: Harmonic band limit [output].
		!!  - N: Azimuthal band limit [output].
		!!  - bl_scoeff: Upper band limit for scaling coefficients [output].
		!!  - alpha: Basis dilation factor [output].
    !!  - filename: Name of the fits input file containing the wavelet
    !!    and scaling coefficients to be read [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - February 2008
    !
    ! Revisions:
    !   February 2008 - Written by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_fits_wav_read_dynamic(wavdyn, scoeff, J, B, N, &
         bl_scoeff, alpha, filename)

			integer, intent(out) :: J
			integer, intent(out) :: B
			integer, intent(out) :: N
			integer, intent(out) :: bl_scoeff
			real(dp), intent(out) :: alpha
      type(ssht_wav_abg), intent(out), allocatable :: wavdyn(:)
			complex(dpc), allocatable, intent(out) :: scoeff(:,:)
      character(len=*), intent(in) :: filename

      character(len=20) :: comment
      integer :: status, unit, blocksize, readwrite
      integer :: nhdu, hdutype, naxis1, naxis2, naxis3
      integer :: hdunum, hdunum_check
      logical :: anynull, file_exists
      integer :: group, dim1, dim2
      real(sp) :: nullval
			integer :: jj, bl_lo, bl_hi, fail
			real(dp), allocatable :: scoeff_temp(:,:)

			! Initialse error status to zero.
      status = 0   

      ! Check if file already exists.
      call ssht_fileio_fits_exists(filename, status, file_exists)
      if(.not. file_exists) then
         call ssht_error(SSHT_ERROR_FILEIO, &
           'ssht_fileio_fits_wav_read', &
           comment_add='File does not exist')
      end if

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Open file as readonly. 
      readwrite = 0    ! Open as readonly.
      call ftopen(unit, filename, readwrite, blocksize, status)

      ! Read size of image from variable size keywords.
      call ftgkyj(unit, 'J', J, comment, status)
      call ftgkyj(unit, 'B', B, comment, status)
      call ftgkyj(unit, 'N', N, comment, status)
      call ftgkyj(unit, 'BLSCOEFF', bl_scoeff, comment, status)
			 call ftgkyd(unit, 'ALPHA', alpha, comment, status)

      ! Check correct number of HDUs in input file.
      ! (First two extension due to primary and dilation, then have one
      ! extension for coefficients at each scale.)
      hdunum = 1 + 2 + (J+1)
      call ftthdu(unit, hdunum_check, status)  ! Number extensions in file.
      if(hdunum_check /= hdunum) then
         call ssht_error(SSHT_ERROR_FILEIO, &
           'ssht_fileio_fits_wav_read', &
           comment_add='Invalid number of headers')
      end if

			! Allocate memory.
      fail = 0
      allocate(wavdyn(0:J), stat=fail)
			allocate(scoeff(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
			allocate(scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1), stat=fail)
			if(fail /= 0) then
				call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, 'ssht_fileio_fits_wav_read')
			end if
			scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = 0d0

			! Read real part of scaling coefficients.
			nhdu = 2 
			call ftmahd(unit, nhdu, hdutype, status)  
			call ftgkyj(unit, 'NAXIS1', naxis1, comment, status)
			call ftgkyj(unit, 'NAXIS2', naxis2, comment, status)
			if(naxis1 /= bl_scoeff .or. naxis2 /= bl_scoeff) then
				call ssht_error(SSHT_ERROR_FILEIO, &
					'ssht_fileio_fits_wav_read', &
					comment_add='Inconsistent scaling coefficients sizes in file')
			end if
			group = 1
			nullval = -999
			dim1 = bl_scoeff
			call ftg2dd(unit, group, nullval, dim1, bl_scoeff, bl_scoeff, &
				scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1), anynull, status)
			scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1)

			! Read imaginary part of scaling coefficients.
			nhdu = 3  
			call ftmahd(unit, nhdu, hdutype, status)  
			call ftgkyj(unit, 'NAXIS1', naxis1, comment, status)
			call ftgkyj(unit, 'NAXIS2', naxis2, comment, status)
			if(naxis1 /= bl_scoeff .or. naxis2 /= bl_scoeff) then
				call ssht_error(SSHT_ERROR_FILEIO, &
					'ssht_fileio_fits_wav_read_dynamic', &
					comment_add='Inconsistent scaling coefficients sizes in file')
			end if
			group = 1
			nullval = -999
			dim1 = bl_scoeff
			call ftg2dd(unit, group, nullval, dim1, bl_scoeff, bl_scoeff, &
				scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1), anynull, status)
			scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) = scoeff(0:bl_scoeff-1, 0:bl_scoeff-1) &
				+ I * scoeff_temp(0:bl_scoeff-1, 0:bl_scoeff-1)

      ! Read wavelet coefficients at each scale.
      do jj = 0,J
	
				! Move to next extension.
				nhdu = nhdu + 1 
				call ftmahd(unit, nhdu, hdutype, status)  

				! Read header and check correct sizes.
				call ftgkyj(unit, 'NAXIS1', naxis1, comment, status)
				call ftgkyj(unit, 'NAXIS2', naxis2, comment, status)
				call ftgkyj(unit, 'NAXIS3', naxis3, comment, status)

        ! Allocate space for wavelet coefficients at given scale.
        allocate(wavdyn(jj)%coeff(0:naxis1-1, 0:naxis2-1, 0:naxis3-1), stat=fail)
        if(fail /= 0) then
           call ssht_error(SSHT_ERROR_MEM_ALLOC_FAIL, &
                'ssht_fileio_fits_wav_read_dynamic')
        end if
        wavdyn(jj)%coeff(0:naxis1-1, 0:naxis2-1, 0:naxis3-1) = 0d0

				! Read coefficients as 3D data cube.
				group = 1
				nullval = -999
				dim1 = naxis1
				dim2 = naxis2
				call ftg3dd(unit, group, nullval, dim1, dim2, &
					naxis1, naxis2, N, &
					wavdyn(jj)%coeff(0:naxis1-1, 0:naxis2-1, 0:naxis3-1), anynull, status)
         if(anynull) then
            call ssht_error(SSHT_ERROR_FILEIO, &
              'ssht_fileio_fits_wav_read_dynamic', &
              comment_add='Null wavelet coefficients read from file')
         end if

      end do

      ! Close fits file.
      call ftclos(unit, status)

      ! Deallocate unit number.
      call ftfiou(unit, status)

      ! Check for errors.
      if (status .gt. 0) call ssht_fileio_fits_error_check(status, .true.)

      ! Free temporary storage used.
      deallocate(scoeff_temp)

    end subroutine ssht_fileio_fits_wav_read_dynamic


    !--------------------------------------------------------------------------
    ! ssht_fileio_fits_error_check
    !
    !! Checks if a fits error has occured and print error message.  Halt
    !! program execution if halt flag is set.
    !!
    !! Variables:
    !!   - status: Fits integer status code [input/output].
    !!   - halt: Logical to indicate whether to halt program execution if an 
    !!     error is detected [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2007
    !
    ! Revisions:
    !   November 2007 - Adapted from CSWT code by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_fits_error_check(status, halt)

      integer, intent(inout) :: status
      logical, intent(in) :: halt

      character(len=30) :: errtext
      character(len=80) :: errmessage

      !  Check if status is OK (no error); if so, simply return.
      if (status .le. 0) return

      ! The FTGERR subroutine returns a descriptive 30-character text 
      ! string that corresponds to the integer error status number.  
      call ftgerr(status,errtext)
      print *,'FITSIO Error Status =',status,': ',errtext

      ! The FTGMSG subroutine retrieves the oldest message from
      ! the stack and shifts any remaining messages on the stack down one
      ! position.  FTGMSG is called repeatedly until a blank message is
      ! returned, which indicates that the stack is empty.  Each error message
      ! may be up to 80 characters in length. 
      call ftgmsg(errmessage)
      do while (errmessage .ne. ' ')
          write(*,*) trim(errmessage)
          call ftgmsg(errmessage)
      end do

      if(halt) stop

    end subroutine ssht_fileio_fits_error_check


    !--------------------------------------------------------------------------
    ! ssht_fileio_fits_exists
    !
    !! Checks if a fits file exists.
    !!
    !! Variables:
    !!   - filename: Name of fits file to check existence of [input].
    !!   - status: Fits integer status code [input].
    !!   - exists: Logical indicating whether the fits file already exists
		!!     [output].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2007
    !
    ! Revisions:
    !   November 2007 - Adapted from CSWT code by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_fits_exists(filename, status, exists)

      character(len=*), intent(in) :: filename
      integer, intent(inout) :: status
      logical, intent(out) :: exists

      integer :: unit, blocksize
      logical :: halt

      ! Simply return if status is already greater than zero.
      if (status .gt. 0) return

      !  Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      call ftopen(unit, filename, 1, blocksize, status)

      ! Check status of opening file.
      if(status == 0) then

        ! File was opened.  Close it and set exists flag accordingly.
        call ftclos(unit, status)
        exists = .true.

      else if (status == 104) then
        
        ! File does not exist.  Reset status and set exists flag accordingly.
         status = 0
         exists = .false.

      else

        ! Some other error occured while opening file.
        halt = .false.
        call ssht_fileio_fits_error_check(status, halt)
        call ftclos(unit, status)
        status = 0
        exists = .true.

      end if

      ! Deallocate unit number.
      call ftfiou(unit, status)

    end subroutine ssht_fileio_fits_exists


    !--------------------------------------------------------------------------
    ! ssht_fileio_fits_del
    !
    !! Deletes a fits file.
    !!
    !! Variables:
    !!   - filename: Name of fits file to detele [input].
    !!   - status: Fits integer status code [input].
    !
    !! @author J. D. McEwen
    !! @version 0.1 - November 2007
    !
    ! Revisions:
    !   November 2007 - Adapted from CSWT code by Jason McEwen
    !--------------------------------------------------------------------------

    subroutine ssht_fileio_fits_del(filename, status)

      character(len=*), intent(in) :: filename
      integer, intent(inout) ::  status

      integer :: unit, blocksize

      ! Simply return if status is greater than zero.
      if (status .gt. 0)return

      ! Get an unused Logical Unit Number to use to open the FITS file.
      call ftgiou(unit,status)

      ! Try to open the file, to see if it exists.
      call ftopen(unit,filename,1,blocksize,status)

      if(status .eq. 0) then
         ! File was opened;  so now delete it.
         call ftdelt(unit,status)
      else if(status .eq. 103) then
         ! File doesn't exist, so just reset status to zero and clear errors.
          status=0
          call ftcmsg
      else
         ! There was some other error opening the file; delete the file anyway.
         status=0
         call ftcmsg
         call ftdelt(unit,status)
      end if

      ! Free the unit number for later reuse.
      call ftfiou(unit, status)

    end subroutine ssht_fileio_fits_del


end module ssht_fileio_mod



