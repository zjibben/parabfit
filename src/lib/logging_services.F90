!!
!! LOGGING_SERVICES
!!
!! This module provides a common set of procedures for writing log/trace
!! messages -- the types of messages typically written to the terminal or
!! simulation log file.  All application code should use the facilities
!! provided here, and eschew use of raw Fortran writes for this purpose.
!!
!! This is simplified version of TRUCHAS_LOGGING_SERVICES.  Parallel-aware
!! aspects have been eliminated, along with features peculiar to Truchas.
!!
!! Neil N. Carlson <nnc@lanl.gov>
!! February 2014
!!
!! PROGRAMMING INTERFACE
!!
!!  There are three categories of log messages. Informational messages are the
!!  typical generic messages, written as-is without any adornment.  These are
!!  logged using LS_INFO.  Warning messages are intended to warn the user of
!!  a possible non-fatal error or exceptional condition.  These are written
!!  prefixed with "Warning:" to be readily distinguished from other messages.
!!  These are logged using LS_WARN. The final category is fatal error messages,
!!  and these are logged using either LS_FATAL or LS_ERROR. The former prefixes
!!  the message with "FATAL:" and gracefully terminates execution. The latter
!!  prefixes the message with "ERROR:" but does not terminate execution.
!!  However in this case it is required that execution ultimately be terminated
!!  using LS_FATAL.  The intent is to permit the logging of multiple fatal
!!  errors before finally terminating execution.
!!
!!  Informational messages are subdivided into several verbosity levels.
!!  Only those messages not exceeding the system verbosity level are logged.
!!  The verbosity levels are provided by the named module constants
!!  LS_VERB_SILENT, LS_VERB_NORMAL, and LS_VERB_NOISY.  The read-only module
!!  variable LS_VERBOSITY gives the system verbosity level.  These values
!!  are of a private derived-type, and thus are only useful for passing to the
!!  procedures described below (the VERBOSITY argument), or in a comparison
!!  expression like "if (LS_VERBOSITY >= LS_VERB_NOISY) ..."  All of the
!!  comparison operators are defined for the type.
!!
!!  What follows is a more detailed description of the procedures.
!!
!!  CALL LS_INITIALIZE (LUNS [,VERBOSITY]) defines the destination for all
!!    log messages and sets the verbosity level.  This must be called before
!!    any of the other subroutines.  LUNS is a rank-1 integer array containing
!!    the logical units where log messages are to be written; they must be
!!    connected for write access.  If not specified, the verbosity level
!!    defaults to LS_VERB_NORMAL.
!!
!!  CALL LS_INFO (MESSAGE [,VERBOSITY]) writes the message passed in the
!!    scalar character argument MESSAGE.  VERBOSITY specifies the verbosity
!!    level of the message, and defaults to LS_VERB_NORMAL.  If the message
!!    verbosity level exceeds the system verbosity level (LS_VERBOSITY) the
!!    message is not written.
!!
!!  CALL LS_INFO (MESSAGE, ADVANCE [,VERBOSITY]) writes the message passed
!!    in the scalar character argument MESSAGE.  If the logical argument
!!    ADVANCE is false, the usual trailing newline is not written, otherwise
!!    the behavior is exactly as for the first variant of LS_INFO.
!!
!!  CALL LS_WARN (MESSAGE) writes the warning message passed in the scalar
!!    character argument MESSAGE.  The written message is prefixed with the
!!    string "Warning: ".
!!
!!  CALL LS_ERROR (MESSAGE) writes the error message passed in the scalar
!!    character argument MESSAGE.  The written message is prefixed with the
!!    string "ERROR: ".  This procedure does not initiate termination of the
!!    simulation, however the caller is expected to ultimately do so using
!!    LS_FATAL.
!!
!!  CALL LS_FATAL (MESSAGE) writes the error message passed in the scalar
!!    character argument MESSAGE.  The written message is prefixed with the
!!    string "FATAL: ".  After logging the message execution is gracefully
!!    terminated with a nonzero exit code.
!!
!!  CALL LS_EXIT logs the fixed normal termination message (no passed message)
!!    and execution is gracefully terminated with a 0 exit code.
!!

#include "f90_assert.fpp"

module logging_services

  implicit none
  private

  public :: LS_initialize
  public :: LS_info, LS_warn, LS_error, LS_fatal, LS_exit

  interface LS_info
    module procedure LS_info_scalar, LS_info_advance
  end interface LS_info

  interface LS_warn
    module procedure LS_warn_scalar
  end interface LS_warn

  interface LS_error
    module procedure LS_error_scalar
  end interface

  !! Log messages are written to these units.
  integer, allocatable, save :: log_unit(:)

  !! Private type to describe verbosity levels.
  type :: verb_level
    private
    integer :: level
  end type verb_level

  !! Named verbosity level constants; normal and noisy will do for now.
  type(verb_level), parameter, public :: LS_VERB_SILENT = verb_level(0)
  type(verb_level), parameter, public :: LS_VERB_NORMAL = verb_level(1)
  type(verb_level), parameter, public :: LS_VERB_NOISY  = verb_level(2)

  !! The system verbosity level; message levels are compared against this.
  type(verb_level), save, public :: LS_verbosity = LS_VERB_NORMAL ! PROTECTED!

  interface operator(==)
    module procedure verb_level_eq
  end interface

  interface operator(/=)
    module procedure verb_level_ne
  end interface

  interface operator(<)
    module procedure verb_level_lt
  end interface

  interface operator(<=)
    module procedure verb_level_le
  end interface

  interface operator(>)
    module procedure verb_level_gt
  end interface

  interface operator(>=)
    module procedure verb_level_ge
  end interface

  public :: operator(==), operator(/=), operator(<), operator(<=), operator(>), operator(>=)

contains

  subroutine LS_initialize (luns, verbosity)
    integer, intent(in) :: luns(:)
    type(verb_level), intent(in), optional :: verbosity
    ASSERT(size(luns) > 0)
    log_unit = luns
    if (present(verbosity)) LS_verbosity = verbosity
  end subroutine LS_initialize

 !!
 !! SPECIFIC LS_INFO PROCEDURES
 !!

  subroutine LS_info_scalar (message, verbosity)
    character(*), intent(in) :: message
    type(verb_level), intent(in), optional :: verbosity
    integer :: n
    type(verb_level) :: level
    if (present(verbosity)) then
      ASSERT(verbosity > LS_VERB_SILENT)
      level = verbosity
    else
      level = LS_VERB_NORMAL
    end if
    if (level <= LS_verbosity) then
      do n = 1, size(log_unit)
        write(log_unit(n),'(a)') message(:len_trim(message))
      end do
    end if
  end subroutine LS_info_scalar

  subroutine LS_info_advance (message, advance, verbosity)
    character(*), intent(in) :: message
    logical, intent(in) :: advance
    type(verb_level), intent(in), optional :: verbosity
    integer :: n
    type(verb_level) :: level
    if (present(verbosity)) then
      ASSERT(verbosity > LS_VERB_SILENT)
      level = verbosity
    else
      level = LS_VERB_NORMAL
    end if
    if (level <= LS_verbosity) then
      if (advance) then
        do n = 1, size(log_unit)
          write(log_unit(n),'(a)') message(:len_trim(message))
        end do
      else
        do n = 1, size(log_unit)
          write(log_unit(n),'(a)',advance='no') message ! intentionally not trimmed
        end do
      end if
    end if
  end subroutine LS_info_advance

 !!
 !! SPECIFIC LS_WARN AND LS_ERROR PROCEDURES
 !!

  subroutine LS_warn_scalar (message)
    character(*), intent(in) :: message
    call labeled_message_scalar ('Warning: ', message)
  end subroutine LS_warn_scalar

  subroutine LS_error_scalar (message)
    character(*), intent(in) :: message
    call labeled_message_scalar ('ERROR: ', message)
  end subroutine LS_error_scalar

  subroutine labeled_message_scalar (label, message)
    character(*), intent(in) :: label, message
    integer :: n
    do n = 1, size(log_unit)
      write(log_unit(n),'(2a)') label, message(:len_trim(message))
    end do
  end subroutine labeled_message_scalar

  subroutine LS_fatal (message)
#ifdef NAGFOR
    use,intrinsic :: f90_unix, only: exit
#endif
    character(*), intent(in) :: message
    character(:), allocatable :: date, time, zone
    call labeled_message_scalar ('FATAL: ', message)
    call timestamp (date, time, zone)
    call LS_info ('Program terminated abnormally on ' // date // ' at ' // time // ' ' // zone)
    call exit (1)
  end subroutine LS_fatal

  subroutine LS_exit
#ifdef NAGFOR
    use,intrinsic :: f90_unix, only: exit
#endif
    character(:), allocatable :: date, time, zone
    call timestamp (date, time, zone)
    call LS_info ('Program terminated normally on ' // date // ' at '// time // ' ' // zone)
    call exit (0)
  end subroutine LS_exit

  subroutine timestamp (date, time, zone)
    character(:), allocatable, intent(out) :: date, time, zone
    character(8)  :: d
    character(10) :: t
    character(5)  :: z
    call date_and_time (date=d, time=t, zone=z)
    date = d(1:4) // '-' // d(5:6) // '-' // d(7:8)
    time = t(1:2) // ':' // t(3:4) // ':' // t(5:6)
    zone = z(1:5)
  end subroutine

 !!
 !! COMPARISON OPERATORS FOR TYPE(VERB_LEVEL) OBJECTS
 !!

  pure logical function verb_level_eq (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_eq = (a%level == b%level)
  end function verb_level_eq

  pure logical function verb_level_ne (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_ne = (a%level /= b%level)
  end function verb_level_ne

  pure logical function verb_level_lt (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_lt = (a%level < b%level)
  end function verb_level_lt

  pure logical function verb_level_le (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_le = (a%level <= b%level)
  end function verb_level_le

  pure logical function verb_level_gt (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_gt = (a%level > b%level)
  end function verb_level_gt

  pure logical function verb_level_ge (a, b)
    type(verb_level), intent(in) :: a, b
    verb_level_ge = (a%level >= b%level)
  end function verb_level_ge

end module logging_services
