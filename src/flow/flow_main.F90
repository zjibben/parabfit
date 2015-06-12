!!
!! FLOW_MAIN
!!
!! A Basic driver for the flow model.
!!
!! Zechariah J. Jibben (zjibben@lanl.gov)
!! June 2015
!!

#include "f90_assert.fpp"

program flow_main

  use,intrinsic :: iso_fortran_env, only: output_unit, error_unit
  use,intrinsic :: iso_c_binding, only: C_NEW_LINE
  use logging_services
  use timer_tree_type
  use parameter_list_type
  use parameter_list_json
  use flow_sim_type
  implicit none

  integer :: n, num_arg, inlun, stat
  character(256) :: arg
  character(:), allocatable :: prog, infile, errmsg
  type(parameter_list), pointer :: params
  type(flow_sim) :: sim

  !! Get the program name from the command line.
  call get_command (arg)
  n = scan(arg, '/', back=.true.)
  prog = trim(arg(n+1:))  ! remove the leading path component, if any

  !! Get the input file name from the command line.
  num_arg = command_argument_count()
  if (num_arg == 1) then
    call get_command_argument (1, arg)
    infile = trim(arg)
  else
    write(error_unit,'(a)') 'usage: ' // prog // ' INFILE'
    stop
  end if

  !! Initialize the logging service routines; output goes to stdout.
  call LS_initialize ([output_unit]) !, LS_VERB_NOISY)

  !! Read the parameter list from the input file.
  open(newunit=inlun,file=infile,action='read',access='stream')
  call parameter_list_from_json_stream (inlun, params, errmsg)
  if (.not.associated(params)) call LS_fatal ("error reading input file:" // C_NEW_LINE // errmsg)
  close(inlun)

  !! Set up the simulation and run it.
  call start_timer ('simulation')
  call sim%init (params)
  call sim%run (stat, errmsg)
  if (stat /= 0) call LS_fatal ('FLOW_SIM: '//errmsg)
  call stop_timer ('simulation')

  !! Write some timing info.
  call LS_info (C_NEW_LINE//'Timing Summary:'//C_NEW_LINE)
  call write_timer_tree (output_unit, indent=3)

  !! And quit.
  call LS_info ('')
  call LS_exit

end program flow_main
