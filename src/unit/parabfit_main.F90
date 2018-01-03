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

  call LS_initialize([output_unit])
  !call set_timer_type (realtime_timing)

  error_cutoff = 1e-5_r8
  do shape = 1,3
    if (shape == 3) error_cutoff = 1e-3_r8
    do mesh = 1,3
      call curvature_grid_refinement_study(mesh, shape, error_cutoff)
    end do
  end do

  call write_timer_tree(output_unit, indent=3)
  call LS_exit()

end program parabfit
