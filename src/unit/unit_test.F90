program unit_test
  use kinds,               only: r8
  use volume_track_module, only: material_flux_unit_test
  use polyhedron_type,     only: polyhedron_unit_test
  use surface_type,        only: surface_unit_test
  use locate_plane_module
  use vof_solver_type,     only: parallel_interfaces_test, intersecting_interfaces_test
  implicit none

  type(locate_plane_hex) :: locate_plane
  real(r8) :: posXflow(6), posxyzn(3), posxyn(3), posxn(3)
  logical  :: test_result

  ! define face velocities [??, ??, -x, +x, ??, ??]
  posXflow = [ 0.0_r8, 0.0_r8, -1.0_r8, 1.0_r8, 0.0_r8, 0.0_r8 ]

  ! define normals
  posxyzn = [1.0_r8, 1.0_r8, 1.0_r8] / sqrt(3.0_r8) ! positive x-y-z
  posxyn  = [1.0_r8, 1.0_r8, 0.0_r8] / sqrt(2.0_r8) ! positive x-y plane
  posxn   = [1.0_r8, 0.0_r8, 0.0_r8]                ! positive x-direction

  ! cube half filled along xyz diagonal
  test_result = locate_plane%unit_test (posxyzn, 0.5_r8, sqrt(3.0_r8)/2.0_r8)
  write(*,*) 'passed locate_plane?     ', test_result

  !test_result = material_flux_unit_test (locate_plane, 0.5_r8*posXflow, 7.0_r8/48.0_r8)

  ! cube half filled along xy diagonal
  ! write(*,*)
  ! test_result = locate_plane%unit_test (posxyn, 0.5_r8, 1.0_r8/sqrt(2.0_r8))
  ! write(*,*) 'passed locate_plane? ', test_result

  !test_result = material_flux_unit_test (locate_plane, 0.5_r8*posXflow, 0.125_r8)
  
  ! cube one eighth filled along xy diagonal
  ! write(*,*)
  ! test_result = locate_plane%unit_test (-posxyn, 0.125_r8, 1.0_r8/(2.0_r8*sqrt(2.0_r8)))
  ! write(*,*) 'passed locate_plane? ', test_result

  !locate_plane%normal = - locate_plane%normal
  ! call material_flux_unit_test (locate_plane,  0.5_r8*posXflow, 0.0_r8)
  ! call material_flux_unit_test (locate_plane, -0.5_r8*posXflow, 0.125_r8)
  
  ! cube 8/10ths filled in x direction
  write(*,*)

  test_result = locate_plane%unit_test (posxn, 0.8_r8, 0.8_r8)
  write(*,*) 'passed locate_plane?     ', test_result
  
  test_result = material_flux_unit_test (locate_plane,  0.25_r8*posXflow, &
       [0.0_r8, 0.0_r8, 0.0_r8, 0.05_r8, 0.0_r8, 0.0_r8])
  write(*,*) 'passed material_flux +x? ', test_result

  test_result = material_flux_unit_test (locate_plane, -0.25_r8*posXflow, &
       [0.0_r8, 0.0_r8, 0.25_r8, 0.0_r8, 0.0_r8, 0.0_r8])
  write(*,*) 'passed material_flux -x? ', test_result


  ! cube 8/10ths filled in -x direction
  write(*,*)
  
  test_result = locate_plane%unit_test (-posxn, 0.8_r8, -0.2_r8)
  write(*,*) 'passed locate_plane?     ', test_result
  
  test_result = material_flux_unit_test (locate_plane,  0.25_r8*posXflow, &
       [0.0_r8, 0.0_r8, 0.0_r8, 0.25_r8, 0.0_r8, 0.0_r8])
  write(*,*) 'passed material_flux +x? ', test_result

  test_result = material_flux_unit_test (locate_plane, -0.25_r8*posXflow, &
       [0.0_r8, 0.0_r8, 0.05_r8, 0.0_r8, 0.0_r8, 0.0_r8])
  write(*,*) 'passed material_flux -x? ', test_result

  ! test polyhedrons and surfaces
  call polyhedron_unit_test ()
  call surface_unit_test ()

  ! multiple interfaces tests
  call parallel_interfaces_test ()
  call intersecting_interfaces_test ()

end program unit_test
