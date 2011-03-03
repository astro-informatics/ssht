! SSHT package to perform spin spherical harmonic transforms
! Copyright (C) 2011  Jason McEwen
! See LICENSE.txt for license details


!------------------------------------------------------------------------------
! ssht_error_mod  -- SSHT library error class
! 
!> Functionality to handle errors that may occur in the ssht library.  Public
!! ssht error codes are defined, with corresponding private error comments and 
!! default halt execution status.
!!
!! @author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
!
! Revisions:
!   October 2010 - Written by Jason McEwen
!------------------------------------------------------------------------------

module ssht_error_mod

  use ssht_types_mod, only: STRING_LEN, SSHT_PROMPT

  implicit none

  private


  !---------------------------------------
  ! Subroutine and function scope
  !---------------------------------------

  public :: ssht_error


  !---------------------------------------
  ! Global variables
  !---------------------------------------

  integer, parameter :: SSHT_ERROR_NUM = 12

  integer, public, parameter :: &
       SSHT_ERROR_NONE = 0, &
       SSHT_ERROR_INIT = 1, &
       SSHT_ERROR_NOT_INIT = 2, &
       SSHT_ERROR_INIT_FAIL = 3, &
       SSHT_ERROR_MEM_ALLOC_FAIL = 4, &
       SSHT_ERROR_ARTH = 5, &
       SSHT_ERROR_SIZE_WARNING = 6, &
       SSHT_ERROR_SIZE_INVALID = 7, &
       SSHT_ERROR_SIZE_NOT_DEF = 8, &
       SSHT_ERROR_ARG_INVALID = 9, &
       SSHT_ERROR_ARG_WARNING = 10, &
       SSHT_ERROR_INDEX_INVALID = 11, &
       SSHT_ERROR_FILEIO = 12

  ! Each element of the error_comment array must have the same length, thus
  ! space with trailing space characters.  When come to use trim to remove 
  ! trailing spaces.
  !> Comment associated with each error type.
  character(len=STRING_LEN), parameter :: &
       error_comment(SSHT_ERROR_NUM+1) = &
       (/ & 
       'No error                                                                 ', &
       'Attempt to initialise object that has already been initialised           ', &
       'Object not initialised                                                   ', &
       'Object initialisation failed                                             ', &
       'Memory allocation failed                                                 ', &
       'Arithmetic exception                                                     ', &
       'Warning: Sizes not in recommended range                                  ', &
       'Invalid sizes                                                            ', &
       'Sizes not defined                                                        ', &
       'Arguments invalid                                                        ', &
       'Argument warning                                                         ', &
       'Index invalid                                                            ', &
       'File IO error                                                            ' &
       /) 

  !> Default program halt status of each error type.
  logical, parameter :: &
       halt_default(SSHT_ERROR_NUM+1) = &
       (/ &
       .false., &
       .true.,  &
       .true.,  &
       .true.,  &
       .true.,  &
       .true.,  &
       .false., &
       .true.,  &
       .true.,  &
       .true.,  &	
       .false.,  &	
       .true.,  &	
       .true.  /)


  !----------------------------------------------------------------------------

contains


  !--------------------------------------------------------------------------
  ! ssht_error
  !
  !> Displays error message corresponding to error_code and halt program 
  !! execution if required.
  !!
  !! Variables:
  !!   - error_code: Integer error code.
  !!   - [procedure]: Procedure name where ssht_error called from.  Displayed 
  !!     when error message printed to screen.
  !!   - [comment_add]: If present, additional comment to append to default 
  !!     error comment.
  !!   - [comment_out]: If present the error comment is copied to comment_out
  !!     on output.
  !!   - [halt_in]: If present overrides default halt value.
  !!
  !! @author <a href="http://www.jasonmcewen.org">Jason McEwen</a>
  !
  ! Revisions:
  !   August 2004 - Written by Jason McEwen
  !--------------------------------------------------------------------------

  subroutine ssht_error(error_code, procedure, comment_add, &
       comment_out, halt_in)

    integer, intent(in) :: error_code
    character(len=*), intent(in), optional :: procedure, comment_add
    character(len=*), intent(inout), optional :: comment_out
    logical, intent(in), optional :: halt_in

    logical :: halt
    character(len=STRING_LEN) :: comment_prefix

    write(comment_prefix, '(a,a)') SSHT_PROMPT, 'SSHT_ERROR: '

    !---------------------------------------
    ! Display error message
    !---------------------------------------

    if(present(procedure)) then

       if(present(comment_add)) then
          write(*,'(a,a,a,a,a,a)') trim(comment_prefix), ' Error ''', &
               trim(error_comment(error_code+1)), &
               ''' occured in procedure ''', &
               trim(procedure), &
               ''''
          write(*,'(a,a,a)') trim(comment_prefix), &
               '  - ', trim(comment_add)
       else
          write(*,'(a,a,a,a,a,a)') trim(comment_prefix), ' Error ''', &
               trim(error_comment(error_code+1)), &
               ''' occured in procedure ''', &
               trim(procedure), &
               ''''
       end if

    else

       if(present(comment_add)) then
          write(*,'(a,a,a)') trim(comment_prefix), &
               ' ', trim(error_comment(error_code+1))
          write(*,'(a,a,a)') trim(comment_prefix), &
               '  - ', trim(comment_add)
       else
          write(*,'(a,a,a)') trim(comment_prefix), ' ', trim(error_comment(error_code+1))
       end if

    end if

    ! Copy error comment if comment_out present.
    if(present(comment_out)) comment_out = error_comment(error_code+1)

    !---------------------------------------
    ! Halt program execution if required
    !---------------------------------------

    if( present(halt_in) ) then
       halt = halt_in
    else
       halt = halt_default(error_code+1)
    end if

    if( halt ) then
       write(*,'(a,a,a,a,a,a)') trim(comment_prefix), ' ', &
            'Halting program execution ', &
            'due to error ''', trim(error_comment(error_code+1)), ''''
       stop
    end if

  end subroutine ssht_error


end module ssht_error_mod
