!!
!! unit_test
!!

program unit_test

  use,intrinsic :: iso_fortran_env, only: output_unit
  use logging_services
  use volume_track_test
  use volume_track_nd_test
  use polygon_type_test
  use polyhedron_type_test
  use multimat_cell_type_test
  use surface_type_test
  use locate_plane_test
  use locate_plane_nd_test
  use vof_solver_type_test
  use analytic_surface_type_test
  use normal_vector_test
  use interface_point_test
  use vof_init_test
  implicit none

  call LS_initialize ([output_unit]) !, LS_VERB_NOISY)

  ! ! vof initialization
  ! call vof_init_test_suite ()

  ! ! locate plane
  ! call locate_plane_unit_test_suite ()

  ! ! fluxing
  ! call volume_track_unit_test_suite ()

  ! ! test polyhedrons and surfaces
  ! call polygon_unit_test ()
  ! call polyhedron_unit_test ()
  ! call surface_unit_test ()

  ! ! nested dissection locate plane
  ! call locate_plane_nd_unit_test_suite ()
  ! call multimat_cell_unit_test_suite ()
  ! call volume_track_nd_unit_test ()

  ! ! multiple interfaces tests
  ! call parallel_interfaces_test ()
  ! call intersecting_interfaces_test ()

  ! analytic surface
  !call analytic_surface_test_suite ()

  call curvature_test_suite ()
  !call normal_vector_test_suite()
  !call interface_point_test_suite()

end program unit_test
