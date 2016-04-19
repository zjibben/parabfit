!! 
!! unit_test
!!

program unit_test

  use volume_track_module,    only: volume_track_unit_test_suite
  use volume_track_nd_module, only: volume_track_nd_unit_test
  use polyhedron_type,        only: polyhedron_unit_test
  use surface_type,           only: surface_unit_test
  use locate_plane_module,    only: locate_plane_unit_test_suite
  use locate_plane_nd_module, only: locate_plane_nd_unit_test_suite
  use vof_solver_type,        only: parallel_interfaces_test, intersecting_interfaces_test
  use multimat_cell_type,     only: multimat_cell_unit_test_suite
  implicit none

  ! locate plane
  call locate_plane_unit_test_suite ()
  
  ! fluxing
  call volume_track_unit_test_suite ()
  
  ! test polyhedrons and surfaces
  call polyhedron_unit_test ()
  call surface_unit_test ()

  ! nested dissection locate plane
  call locate_plane_nd_unit_test_suite ()
  call multimat_cell_unit_test_suite ()
  call volume_track_nd_unit_test ()

  ! multiple interfaces tests
  call parallel_interfaces_test ()
  call intersecting_interfaces_test ()

end program unit_test
