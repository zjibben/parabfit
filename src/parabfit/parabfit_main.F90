program parabfit

  use, intrinsic :: iso_fortran_env, only: output_unit
  use kinds, only: r8
  use logging_services
  use unstr_mesh_type
  use mesh_geom_type
  use timer_tree_type
  use curvature_tests
  implicit none

  real(r8) :: error_cutoff
  integer :: mesh, shape
  logical :: run_full
  character(16) :: arg

  call LS_initialize([output_unit])
  !call set_timer_type (realtime_timing)

  ! check if user requested full run
  call get_command_argument(1, arg)
  run_full = trim(adjustl(arg))=='full'

  ! make required dump directories if they don't already exist
  call execute_command_line('mkdir -p conv')
  call execute_command_line('mkdir -p fields')

  ! run
  if (run_full) then
    error_cutoff = 1e-5_r8
    do shape = 1,3
      if (shape == 3) error_cutoff = 1e-3_r8
      do mesh = 1,3
        call curvature_grid_refinement_study(mesh, shape, error_cutoff)
      end do
    end do
  else
    call ellipsoid_coarse()
  end if

  print *
  call write_timer_tree(output_unit, indent=3)
  call LS_exit()

end program parabfit

subroutine ellipsoid_coarse()

  use kinds, only: r8
  use curvature_tests

  type(curvature_test) :: test_params
  real(r8) :: lnormFT(3), lnormHF(3)

  call test_params%init(TET, ELLIPSOID, 1e-5_r8)
  call mesh_test(test_params, 20, lnormFT, lnormHF, 0, 0)

end subroutine ellipsoid_coarse
